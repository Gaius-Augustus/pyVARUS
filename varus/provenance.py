"""What is needed to rebuild VARUS.bam without keeping it.

``varus run`` writes two small files next to ``VARUS.bam``:

``VARUS.manifest.tsv``
    One row per batch in the BAM: SRA run, spot range, paired flag, aligner
    preset, the version of the splice-site DB the batch was aligned with and
    the aligner's thread count (HISAT2 places some reads differently with a
    different ``-p``). The header records the varus and tool versions, the genome's MD5
    and the command line.

``VARUS.splicedb.log.gz``
    Every change of the aligner's splice-site DB (``intronDB.splice_sites``
    or ``intronDB.junc.bed``) as ``version<TAB>+|-<TAB>line``. The DB grows
    while VARUS runs and each batch is aligned against the DB of that
    moment, so the final ``introns.gff`` cannot stand in for it.

Together with the genome FASTA they are all ``varus replay``
(:mod:`varus.replay`) needs to download the same reads again and align them
the same way.
"""

from __future__ import annotations

import datetime as _dt
import gzip
import hashlib
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

MANIFEST_NAME = "VARUS.manifest.tsv"
SPLICE_LOG_NAME = "VARUS.splicedb.log.gz"
FORMAT_VERSION = "1"

MANIFEST_COLUMNS = [
    "batch", "accession", "n", "x", "paired", "platform", "preset",
    "n_batches", "db_version", "align_threads", "uniq_pct",
]


# ---------------------------------------------------------------------------
# Checksums and tool versions
# ---------------------------------------------------------------------------

def file_md5(path: Path, chunk: int = 1 << 22) -> str:
    """MD5 of a file, or "" if it cannot be read."""
    h = hashlib.md5()
    try:
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(chunk), b""):
                h.update(block)
    except OSError:
        return ""
    return h.hexdigest()


def tool_version(exe: str) -> str:
    """First non-empty line of ``<exe> --version``, or "not found"."""
    path = shutil.which(exe)
    if path is None:
        return "not found"
    try:
        p = subprocess.run([path, "--version"], capture_output=True, text=True,
                           errors="replace", timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    for line in (p.stdout + "\n" + p.stderr).splitlines():
        if line.strip():
            return line.strip()
    return "unknown"


# ---------------------------------------------------------------------------
# Splice-DB change log
# ---------------------------------------------------------------------------

def _db_line_key(line: str) -> tuple:
    """Sort key of a splice-DB line: chrom, then numeric start and end.

    Both writers (:func:`varus.strand.write_hisat2_splice_sites` and
    :func:`varus.strand.write_minimap2_junc_bed`) emit lines in this order.
    """
    f = line.split("\t")
    try:
        return (f[0], int(f[1]), int(f[2]), f[3:])
    except (IndexError, ValueError):
        return (line, 0, 0, [])


class SpliceDBLog:
    """Append-only record of every version of the aligner's splice-site DB.

    :meth:`record` is called after each (re)write of the DB file; it diffs the
    file against the previous version and appends the added and removed
    lines under a new version number. Version 0 is "no DB". Each version is
    appended as its own gzip member, so the log stays readable if the run
    dies.
    """

    def __init__(self, path: Path, fmt: str) -> None:
        self.path = Path(path)
        self.version = 0
        self._lines: Set[str] = set()
        with gzip.open(self.path, "wt", encoding="utf-8") as f:
            f.write(f"#varus_splice_db_log={FORMAT_VERSION}\n#format={fmt}\n")

    def record(self, db_path: Path) -> int:
        """Log the current content of ``db_path``; return its version."""
        try:
            with open(db_path, encoding="utf-8") as f:
                new = {ln.rstrip("\n") for ln in f if ln.strip()}
        except OSError:
            new = set()
        added = new - self._lines
        removed = self._lines - new
        if not added and not removed:
            return self.version
        self.version += 1
        v = self.version
        with gzip.open(self.path, "at", encoding="utf-8") as f:
            for ln in sorted(removed, key=_db_line_key):
                f.write(f"{v}\t-\t{ln}\n")
            for ln in sorted(added, key=_db_line_key):
                f.write(f"{v}\t+\t{ln}\n")
        self._lines = new
        return v


class SpliceDBSnapshots:
    """Rebuild the splice-site DB of any version from a :class:`SpliceDBLog`."""

    def __init__(self, path: Path) -> None:
        self.fmt = ""
        self._ops: List[Tuple[int, str, str]] = []
        with gzip.open(path, "rt", encoding="utf-8") as f:
            for raw in f:
                raw = raw.rstrip("\n")
                if raw.startswith("#"):
                    if raw.startswith("#format="):
                        self.fmt = raw.split("=", 1)[1]
                    continue
                if not raw:
                    continue
                v, op, line = raw.split("\t", 2)
                self._ops.append((int(v), op, line))
        self._ops.sort(key=lambda t: t[0])     # stable: keeps - before + per version
        self.max_version = self._ops[-1][0] if self._ops else 0
        self._lines: Set[str] = set()
        self._pos = 0
        self._version = 0

    def lines(self, version: int) -> List[str]:
        """The DB lines of ``version``, in the order the writers emit them."""
        if version > self.max_version:
            raise ValueError(f"splice DB version {version} not in the log "
                             f"(last version {self.max_version})")
        if version < self._version:
            self._lines, self._pos, self._version = set(), 0, 0
        while self._pos < len(self._ops) and self._ops[self._pos][0] <= version:
            _, op, line = self._ops[self._pos]
            if op == "+":
                self._lines.add(line)
            else:
                self._lines.discard(line)
            self._pos += 1
        self._version = version
        return sorted(self._lines, key=_db_line_key)

    def write(self, version: int, out: Path) -> Optional[Path]:
        """Write the DB of ``version`` to ``out``; ``None`` for version 0 (no DB)."""
        if version <= 0:
            return None
        lines = self.lines(version)
        with open(out, "w", encoding="utf-8") as f:
            for ln in lines:
                f.write(ln + "\n")
        return out


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def write_manifest(path: Path, header: Dict[str, str], rows: List[dict]) -> Path:
    """Write ``VARUS.manifest.tsv`` (replaced atomically)."""
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("# VARUS manifest: the read batches in VARUS.bam and how each was "
                "aligned.\n")
        f.write("# Keep this file and " + SPLICE_LOG_NAME + " (and the genome) to "
                "rebuild the BAM with `varus replay`.\n")
        for k, v in header.items():
            v = str(v).replace("\n", " ").replace("\t", " ")
            f.write(f"#{k}={v}\n")
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in MANIFEST_COLUMNS) + "\n")
    os.replace(tmp, path)
    return path


def read_manifest(path: Path) -> Tuple[Dict[str, str], List[dict]]:
    """Parse ``VARUS.manifest.tsv`` into (header, rows)."""
    header: Dict[str, str] = {}
    rows: List[dict] = []
    cols: Optional[List[str]] = None
    with open(path, encoding="utf-8") as f:
        for raw in f:
            raw = raw.rstrip("\n")
            if raw.startswith("# ") or not raw.strip():
                continue
            if raw.startswith("#"):
                k, _, v = raw[1:].partition("=")
                header[k] = v
                continue
            if cols is None:
                cols = raw.split("\t")
                missing = set(MANIFEST_COLUMNS) - set(cols)
                if missing:
                    raise ValueError(f"{path}: missing columns {sorted(missing)}")
                continue
            d = dict(zip(cols, raw.split("\t")))
            rows.append({
                "batch": int(d["batch"]),
                "accession": d["accession"],
                "n": int(d["n"]),
                "x": int(d["x"]),
                "paired": d["paired"] == "1",
                "platform": d["platform"],
                "preset": d["preset"],
                "n_batches": int(d["n_batches"]),
                "db_version": int(d["db_version"]),
                "align_threads": int(d["align_threads"]),
                "uniq_pct": float(d["uniq_pct"]),
            })
    if cols is None:
        raise ValueError(f"{path}: no column header; not a VARUS manifest")
    return header, rows


def now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
