"""Sample batches from local gzipped FASTQ files instead of (or next to) SRA.

A local library is handed to ``varus run --fastq R1.fq.gz[,R2.fq.gz]`` and
takes part in the online sampling loop like an SRA run: it has a number of
spots, is cut into batches of ``batch_size`` spots, and the controller
"downloads" batch ``k`` by extracting spots ``[k*bs, (k+1)*bs)`` as FASTA
into the usual ``batches/<name>/N<n>X<x>/`` directory.

Random access into a gzip stream
--------------------------------
gzip has no random access, and the loop draws the batches of a run in a
shuffled order (a FASTQ derived from a coordinate-sorted BAM would otherwise
be sampled from one end of the genome only). Rather than decompressing the
file once per batch or writing a re-encoded copy to disk, the file is read
*once* when the run starts (:class:`GzipFastqIndex`); that pass

* counts the records and bases (``total_spots``, ``avg_len``),
* remembers the uncompressed byte offset of every ``stride``-th record, and
* keeps a copy of the zlib decompressor state every ``snapshot_bytes`` of
  uncompressed data (``zlib.decompressobj().copy()``: ~40 KB each).

Extracting a batch then seeks to the nearest snapshot, resumes the
decompressor from its copy, skips to the record offset and parses the
records. One batch costs at most ``snapshot_bytes`` of extra decompression
(~0.1 s). The index lives in memory only; ``varus replay`` rebuilds it.

Input files must be gzip-compressed FASTQ (``.fastq.gz`` or ``.fq.gz``;
multi-member gzip as written by ``pigz``, ``bgzip`` or Illumina instruments is
fine). Paired-end libraries are given as two files with the mates in the same
order; interleaved files are not supported.
"""

from __future__ import annotations

import bisect
import logging
import os
import re
import zlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple

import numpy as np

from varus.download import BatchPaths, batch_dir_for
from varus.runlist import RunRecord

log = logging.getLogger(__name__)

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")
_GZIP_MAGIC = b"\x1f\x8b"

# Compressed bytes fed to the decompressor per step.
_CHUNK = 1 << 20
# Uncompressed bytes between two decompressor snapshots (extra decompression
# per extracted batch is at most this much).
DEFAULT_SNAPSHOT_BYTES = 32 << 20
# Every ``stride``-th record start is remembered (8 bytes each).
DEFAULT_STRIDE = 1024


class LocalReadsError(RuntimeError):
    """A local FASTQ file is unusable (not gzip, not FASTQ, mates differ)."""


# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------

@dataclass
class _Snapshot:
    comp_offset: int            # file offset of the next compressed chunk
    uncomp_offset: int          # uncompressed bytes produced so far
    state: Optional[object]     # zlib decompressobj copy; None = fresh stream


class GzipFastqIndex:
    """Random access to the records of one gzipped FASTQ file (see module doc)."""

    def __init__(self, path: Path, *, stride: int = DEFAULT_STRIDE,
                 snapshot_bytes: int = DEFAULT_SNAPSHOT_BYTES,
                 chunk_bytes: int = _CHUNK) -> None:
        self.path = Path(path)
        self.stride = max(1, int(stride))
        self.chunk = max(1024, int(chunk_bytes))
        self.snapshot_bytes = max(self.chunk, int(snapshot_bytes))
        self.n_records = 0
        self.n_bases = 0
        self.uncompressed_bytes = 0
        self._starts: np.ndarray = np.zeros(0, dtype=np.int64)
        self._snapshots: List[_Snapshot] = []
        self._build()

    # -- building ----------------------------------------------------------

    def _build(self) -> None:
        check_fastq_gz(self.path)
        starts: List[int] = []
        snaps: List[_Snapshot] = [_Snapshot(0, 0, None)]
        d = zlib.decompressobj(wbits=31)
        n_lines = 0           # complete lines seen so far
        last_nl = -1          # uncompressed offset of the last newline seen
        uncomp = 0            # uncompressed bytes produced so far
        comp = 0              # compressed bytes consumed so far
        bases = 0
        since_snap = 0
        last_byte = 0         # last uncompressed byte of the previous chunk
        stride = self.stride
        with open(self.path, "rb") as fh:
            while True:
                chunk = fh.read(self.chunk)
                if not chunk:
                    break
                comp += len(chunk)
                out, d = _inflate(d, chunk, self.path)
                if out:
                    base = uncomp
                    arr = np.frombuffer(out, dtype=np.uint8)
                    nl = np.flatnonzero(arr == 10)
                    if nl.size:
                        pos = nl.astype(np.int64) + base
                        idx = np.arange(n_lines, n_lines + nl.size, dtype=np.int64)
                        prev = np.empty_like(pos)
                        prev[0] = last_nl
                        prev[1:] = pos[:-1]
                        # sequence lines (2nd of each record): sum their lengths,
                        # not counting a '\r' before the newline (CRLF files)
                        seq = (idx & 3) == 1
                        bases += int((pos[seq] - prev[seq] - 1).sum())
                        nl_seq = nl[seq]
                        before = np.where(nl_seq > 0, arr[np.maximum(nl_seq - 1, 0)],
                                          last_byte)
                        bases -= int(np.count_nonzero(before == 13))
                        # a record starts after the newline closing a 4th line
                        rec_end = (idx & 3) == 3
                        rec_idx = (idx[rec_end] + 1) >> 2        # index of the next record
                        keep = (rec_idx % stride) == 0
                        starts.extend((pos[rec_end][keep] + 1).tolist())
                        n_lines += int(nl.size)
                        last_nl = int(pos[-1])
                    uncomp += len(out)
                    since_snap += len(out)
                    last_byte = out[-1]
                if since_snap >= self.snapshot_bytes:
                    snaps.append(_Snapshot(comp, uncomp, d.copy() if d is not None else None))
                    since_snap = 0
        if d is not None:
            raise LocalReadsError(
                f"{self.path}: gzip stream is truncated or corrupt (no end-of-stream "
                "marker)")
        # an unterminated last line still counts
        if uncomp > 0 and last_nl != uncomp - 1:
            n_lines += 1
            if (n_lines - 1) & 3 == 1:
                bases += uncomp - last_nl - 1
        if n_lines % 4:
            raise LocalReadsError(
                f"{self.path}: {n_lines} lines is not a multiple of 4; not a "
                "4-line-per-record FASTQ file")
        self.n_records = n_lines // 4
        self.n_bases = bases
        self.uncompressed_bytes = uncomp
        if self.n_records == 0:
            raise LocalReadsError(f"{self.path}: no reads")
        # record 0 starts at 0; the file may hold more starts than records
        # (a start after the final record is the EOF offset)
        arr = np.array([0] + starts, dtype=np.int64)
        n_keep = (self.n_records - 1) // stride + 1
        self._starts = arr[:n_keep]
        self._snapshots = snaps

    # -- extraction --------------------------------------------------------

    @property
    def avg_len(self) -> float:
        return self.n_bases / self.n_records if self.n_records else 0.0

    def records(self, n: int, x: int) -> Iterator[Tuple[bytes, bytes]]:
        """Yield ``(name, sequence)`` of spots ``n..x`` (inclusive, 0-based)."""
        if n < 0 or n >= self.n_records:
            raise ValueError(f"spot {n} outside 0..{self.n_records - 1}")
        x = min(x, self.n_records - 1)
        want = x - n + 1
        if want <= 0:
            return
        k = n // self.stride
        offset = int(self._starts[k])
        skip = n - k * self.stride           # records to skip after `offset`
        stream = self._stream_from(offset)
        try:
            yield from self._parse(stream, n, skip, want)
        finally:
            stream.close()

    def _parse(self, stream: Iterator[bytes], n: int, skip: int, want: int
               ) -> Iterator[Tuple[bytes, bytes]]:
        buf = bytearray()
        line = 0                              # line within the current record
        got = 0
        name = seq = b""
        eof = False
        while got < want:
            # ensure a complete line is in buf
            nl = buf.find(b"\n")
            while nl < 0 and not eof:
                piece = next(stream, None)
                if piece is None:
                    eof = True
                    break
                buf += piece
                nl = buf.find(b"\n")
            if nl < 0:
                if not buf:
                    raise LocalReadsError(
                        f"{self.path}: ended after {n + got} of {self.n_records} reads")
                text, buf = bytes(buf), bytearray()   # unterminated last line
            else:
                text = bytes(buf[:nl])
                del buf[:nl + 1]
            if text.endswith(b"\r"):
                text = text[:-1]
            if line == 0:
                if skip == 0:
                    if not text.startswith(b"@"):
                        raise LocalReadsError(
                            f"{self.path}: read {n + got} does not start with '@'")
                    name = text[1:].split(None, 1)[0] if len(text) > 1 else b""
            elif line == 1:
                seq = text
            elif line == 2:
                if skip == 0 and not text.startswith(b"+"):
                    raise LocalReadsError(
                        f"{self.path}: read {n + got}: third line does not start with '+'")
            else:
                if skip > 0:
                    skip -= 1
                else:
                    if len(text) != len(seq):
                        raise LocalReadsError(
                            f"{self.path}: read {n + got}: quality and sequence "
                            "lengths differ (truncated or not FASTQ)")
                    yield name, seq
                    got += 1
            line = (line + 1) & 3

    def _stream_from(self, offset: int) -> Iterator[bytes]:
        """Uncompressed bytes from ``offset`` on (resumed at the nearest snapshot)."""
        i = bisect.bisect_right([s.uncomp_offset for s in self._snapshots], offset) - 1
        snap = self._snapshots[max(0, i)]
        d = snap.state.copy() if snap.state is not None else zlib.decompressobj(wbits=31)
        to_skip = offset - snap.uncomp_offset
        with open(self.path, "rb") as fh:
            fh.seek(snap.comp_offset)
            while True:
                chunk = fh.read(self.chunk)
                if not chunk:
                    return
                out, d = _inflate(d, chunk, self.path)
                if to_skip:
                    if len(out) <= to_skip:
                        to_skip -= len(out)
                        continue
                    out, to_skip = out[to_skip:], 0
                if out:
                    yield out

    def write_fasta(self, n: int, x: int, out: Path) -> int:
        """Write spots ``n..x`` as FASTA to ``out``; return the number written."""
        tmp = out.with_name(out.name + ".tmp")
        count = 0
        with open(tmp, "wb") as f:
            for name, seq in self.records(n, x):
                count += 1
                f.write(b">" + (name or str(n + count - 1).encode()) + b"\n" + seq + b"\n")
        os.replace(tmp, out)
        return count


def _inflate(d, chunk: bytes, path: Path):
    """Feed ``chunk`` to decompressor ``d``; handle gzip member boundaries.

    Returns ``(output, decompressor)``; the decompressor is ``None`` when a
    gzip member ended exactly at the end of ``chunk`` (a fresh one is started
    on the next chunk).
    """
    out = bytearray()
    data = chunk
    try:
        while True:
            if d is None:
                d = zlib.decompressobj(wbits=31)
            out += d.decompress(data)
            if not d.eof:
                return bytes(out), d
            rest = d.unused_data
            d = None
            if not rest:
                return bytes(out), None
            # another member follows (pigz, bgzip, Illumina): skip padding
            rest = rest.lstrip(b"\x00")
            if not rest:
                return bytes(out), None
            data = rest
    except zlib.error as e:
        raise LocalReadsError(f"{path}: gzip stream is corrupt: {e}") from e


# ---------------------------------------------------------------------------
# Validation and naming
# ---------------------------------------------------------------------------

def check_fastq_gz(path: Path) -> None:
    """Raise :class:`LocalReadsError` unless ``path`` is a gzipped FASTQ file."""
    path = Path(path)
    name = path.name.lower()
    if not name.endswith(FASTQ_SUFFIXES):
        raise LocalReadsError(
            f"{path}: local read files must be gzipped FASTQ named *.fastq.gz "
            "or *.fq.gz")
    if not path.is_file():
        raise LocalReadsError(f"{path}: file not found")
    with open(path, "rb") as fh:
        head = fh.read(2)
    if head != _GZIP_MAGIC:
        raise LocalReadsError(f"{path}: not gzip-compressed (expected a .fastq.gz / "
                              ".fq.gz file; compress it with gzip or pigz)")


def strip_fastq_suffix(name: str) -> str:
    low = name.lower()
    for suf in FASTQ_SUFFIXES:
        if low.endswith(suf):
            return name[:-len(suf)]
    return name


_MATE_SUFFIX = re.compile(r"[._-]?(R?[12]|R?[12]_001)$", re.IGNORECASE)
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


def run_name(r1: Path, r2: Optional[Path] = None) -> str:
    """Name of a local run: the file name without suffix and mate marker.

    ``lib_1.fq.gz`` + ``lib_2.fq.gz`` -> ``lib``; ``S1_R1_001.fastq.gz`` +
    ``S1_R2_001.fastq.gz`` -> ``S1``; a single-end ``reads.fq.gz`` -> ``reads``.
    Characters outside ``[A-Za-z0-9._-]`` become ``_`` (the name is a
    directory name under ``batches/`` and a column in the manifest).
    """
    s1 = strip_fastq_suffix(Path(r1).name)
    if r2 is None:
        name = s1
    else:
        s2 = strip_fastq_suffix(Path(r2).name)
        m1, m2 = _MATE_SUFFIX.search(s1), _MATE_SUFFIX.search(s2)
        if m1 and m2 and s1[:m1.start()] == s2[:m2.start()] and m1.start() > 0:
            name = s1[:m1.start()]
        else:
            # longest common prefix, trimmed of separators
            i = 0
            while i < min(len(s1), len(s2)) and s1[i] == s2[i]:
                i += 1
            name = s1[:i].rstrip("._-") or s1
    name = _UNSAFE.sub("_", name).strip("._") or "local"
    return name


def parse_fastq_specs(specs: Sequence[str]) -> List[Tuple[Path, Optional[Path]]]:
    """``["a.fq.gz", "b_1.fq.gz,b_2.fq.gz"]`` -> ``[(a, None), (b_1, b_2)]``."""
    out: List[Tuple[Path, Optional[Path]]] = []
    for spec in specs:
        parts = [p for p in spec.split(",") if p.strip()]
        if not parts or len(parts) > 2:
            raise LocalReadsError(
                f"--fastq {spec!r}: give one file (single-end) or two files "
                "separated by a comma (paired-end)")
        r1 = Path(parts[0]).expanduser()
        r2 = Path(parts[1]).expanduser() if len(parts) == 2 else None
        if r2 is not None and r1.resolve() == r2.resolve():
            raise LocalReadsError(f"--fastq {spec!r}: the two mate files are the same file")
        out.append((r1, r2))
    return out


# ---------------------------------------------------------------------------
# Local runs
# ---------------------------------------------------------------------------

@dataclass
class LocalRun:
    """One local library: its files and their indexes."""
    name: str
    r1: Path
    r2: Optional[Path]
    idx1: GzipFastqIndex
    idx2: Optional[GzipFastqIndex] = None

    @property
    def paired(self) -> bool:
        return self.r2 is not None

    @property
    def total_spots(self) -> int:
        return self.idx1.n_records

    @property
    def avg_len(self) -> float:
        """Average spot length (both mates for paired-end, like SRA's avg_len)."""
        bases = self.idx1.n_bases + (self.idx2.n_bases if self.idx2 else 0)
        return bases / self.total_spots if self.total_spots else 0.0

    @property
    def total_bases(self) -> int:
        return self.idx1.n_bases + (self.idx2.n_bases if self.idx2 else 0)

    def files(self) -> List[Path]:
        return [self.r1] if self.r2 is None else [self.r1, self.r2]

    def record(self, platform: str = "") -> RunRecord:
        return RunRecord(
            accession=self.name, total_spots=self.total_spots,
            total_bases=self.total_bases, avg_len=self.avg_len,
            paired=self.paired, colorspace=False, platform=platform,
            bioproject="local",
        )


def index_local_runs(specs: Sequence[Tuple[Path, Optional[Path]]], *,
                     threads: int = 4, names: Optional[Sequence[str]] = None,
                     **index_kw) -> List[LocalRun]:
    """Index all libraries (one thread per file); names must be unique."""
    files = [(Path(r1), Path(r2) if r2 is not None else None) for r1, r2 in specs]
    names = list(names) if names is not None else [run_name(r1, r2) for r1, r2 in files]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise LocalReadsError(
            f"--fastq: several libraries would get the same name {sorted(dup)}; "
            "rename the files so their names (without .fq.gz and the mate "
            "marker _1/_2) differ")
    paths: List[Path] = []
    for r1, r2 in files:
        paths += [r1] if r2 is None else [r1, r2]
    for p in paths:
        check_fastq_gz(p)
    log.info("Indexing %d local FASTQ file(s) for %d run(s); this reads each file once",
             len(paths), len(files))
    workers = max(1, min(int(threads), len(paths)))
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="varus-fqidx") as ex:
        futs = {p: ex.submit(GzipFastqIndex, p, **index_kw) for p in paths}
        idx = {p: f.result() for p, f in futs.items()}
    runs: List[LocalRun] = []
    for (r1, r2), name in zip(files, names):
        idx1 = idx[r1]
        idx2 = idx[r2] if r2 is not None else None
        if idx2 is not None and idx1.n_records != idx2.n_records:
            raise LocalReadsError(
                f"{name}: mate files hold different numbers of reads "
                f"({r1}: {idx1.n_records}, {r2}: {idx2.n_records})")
        run = LocalRun(name=name, r1=r1, r2=r2, idx1=idx1, idx2=idx2)
        log.info("Local run %s: %d spots, avg length %.1f, %s", name, run.total_spots,
                 run.avg_len, "paired" if run.paired else "single-end")
        runs.append(run)
    return runs


def extract_local_batch(run: LocalRun, n: int, x: int, outdir: Path) -> BatchPaths:
    """Write spots ``n..x`` of ``run`` as FASTA into the batch directory.

    The layout and file names match :func:`varus.download.download_batch`
    (``<acc>.fasta`` or ``<acc>_1.fasta`` + ``<acc>_2.fasta``), so the rest of
    the controller does not care where a batch came from.
    """
    bdir = batch_dir_for(outdir, run.name, n, x)
    bdir.mkdir(parents=True, exist_ok=True)
    try:
        if run.idx2 is None:
            r1 = bdir / f"{run.name}.fasta"
            run.idx1.write_fasta(n, x, r1)
            return BatchPaths(r1=r1, r2=None, batch_dir=bdir)
        r1 = bdir / f"{run.name}_1.fasta"
        r2 = bdir / f"{run.name}_2.fasta"
        c1 = run.idx1.write_fasta(n, x, r1)
        c2 = run.idx2.write_fasta(n, x, r2)
        if c1 != c2:
            raise LocalReadsError(f"{run.name}: mates {n}-{x} differ in read count "
                                  f"({c1} vs {c2})")
        return BatchPaths(r1=r1, r2=r2, batch_dir=bdir)
    except (OSError, ValueError, LocalReadsError) as e:
        raise RuntimeError(f"extracting spots {n}-{x} of {run.name} failed: {e}") from e
