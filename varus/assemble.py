"""StringTie assembly and stranded intron hints from ``VARUS.bam``.

With these two files, the BAM can be deleted: they are everything
Paludamentum reads from a VARUS BAM (see ``docs/plan_stringtie_assembly.md``).

``stringtie.gtf``
    ``stringtie -p N -o stringtie.gtf VARUS.bam`` (``-L`` for long reads;
    ``--mix short.bam long.bam`` for both), without other options. That is
    how Paludamentum assembles a BAM.

``hints.gff``
    What ``bam2hints --intronsonly`` followed by
    ``filterIntronsFindStrand.pl genome --score`` writes for the same BAM.
    :func:`bam2hints_introns` ports bam2hints' intron rules (AUGUSTUS,
    ``auxprogs/bam2hints/bam2hints.cc``). One difference on purpose:
    bam2hints keeps the multiplicity in an ``unsigned short`` and wraps
    above 65535 (65537 identical introns come out as ``mult`` 1); here the
    count is exact.

    :mod:`varus.introns` is a different, simpler extraction (every CIGAR
    ``N``, no length window) that feeds the aligner's splice-site DB and
    ``introns.gff``; it is not changed.
"""

from __future__ import annotations

import hashlib
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from varus.provenance import tool_version

log = logging.getLogger(__name__)

STRINGTIE_NAME = "stringtie.gtf"
# Exit status of `varus run` / `varus replay` when StringTie or the hint
# extraction failed; VARUS.bam is kept (also with --drop-bam).
EXIT_ASSEMBLY_FAILED = 4
HINTS_NAME = "hints.gff"

# bam2hints defaults (all left at their defaults by Paludamentum)
MAX_GAP_LEN = 14          # gaps up to this length are closed (deletions)
MIN_INTRON_LEN = 32       # gaps between MAX_GAP_LEN and this drop the alignment
MAX_INTRON_LEN = 350_000  # longer gaps drop the alignment
MAX_QGAP_LEN = 5          # query gap after a block; longer: no intron there
MIN_END_BLOCK_LEN = 8     # minimal length of the first and last block
MAX_GENE_LEN = 400_000    # alignments spanning more are ignored
PRIORITY = 4
SOURCE = "E"

_MATCH = (0, 7, 8)        # M, =, X
_REF_GAP = (2, 3)         # D, N
_IGNORED = (4, 5, 6)      # S, H, P
_QUERY_GAP = 1            # I


def require_stringtie(stringtie: str = "stringtie") -> None:
    """Fail before a long run starts, not at its end."""
    if shutil.which(stringtie) is None:
        raise SystemExit(
            f"{stringtie} not found on PATH. pyVARUS assembles VARUS.bam with "
            "StringTie 3.0.3 (the version Paludamentum uses); install it from "
            "https://github.com/gpertea/stringtie/releases or use the pyVARUS container."
        )


# ---------------------------------------------------------------------------
# Intron hints (bam2hints --intronsonly)
# ---------------------------------------------------------------------------

def _alignment_introns(start0: int, cigar) -> List[Tuple[int, int]]:
    """Intron hints of one alignment, as bam2hints derives them.

    ``start0`` is the 0-based alignment start, ``cigar`` pysam's
    ``cigartuples``. Returns 1-based inclusive ``(start, end)`` pairs.
    """
    # PSL-like blocks: length, 1-based query start, 1-based target start
    b: List[int] = []
    q: List[int] = []
    t: List[int] = []
    qoff, toff = 1, start0 + 1
    for op, ln in cigar:
        if op in _MATCH:
            if b and t[-1] + b[-1] == toff and q[-1] + b[-1] == qoff:
                b[-1] += ln
            else:
                b.append(ln)
                q.append(qoff)
                t.append(toff)
            qoff += ln
            toff += ln
        elif op in _REF_GAP:
            toff += ln
        elif op == _QUERY_GAP:
            qoff += ln
        elif op in _IGNORED:
            pass
        else:                       # B or unknown: bam2hints drops the alignment
            return []
    n = len(b)
    if n < 2 or t[-1] + b[-1] - t[0] > MAX_GENE_LEN:
        return []

    # filter blocks as blat2hints.pl: close short gaps, keep intron-sized ones
    begins: List[int] = []
    ends: List[int] = []
    fol_ok: List[bool] = []
    for i in range(n):
        gap = MIN_INTRON_LEN if not begins else t[i] - ends[-1] - 1
        ok = i < n - 1 and q[i + 1] - q[i] - b[i] <= MAX_QGAP_LEN
        if MIN_INTRON_LEN <= gap <= MAX_INTRON_LEN:
            begins.append(t[i])
            ends.append(t[i] + b[i] - 1)
            fol_ok.append(ok)
        elif gap <= MAX_GAP_LEN:
            ends[-1] = t[i] + b[i] - 1
            fol_ok[-1] = ok
        else:
            return []

    m = len(begins)
    out: List[Tuple[int, int]] = []
    for i in range(m - 1):
        if not fol_ok[i]:
            continue
        if i == 0 and ends[0] - begins[0] + 1 < MIN_END_BLOCK_LEN:
            continue
        if i == m - 2 and ends[i + 1] - begins[i + 1] + 1 < MIN_END_BLOCK_LEN:
            continue
        out.append((ends[i] + 1, begins[i + 1] - 1))
    return out


def bam2hints_introns(bam_path: Path, threads: int = 1) -> Dict[str, Dict[Tuple[int, int], int]]:
    """``bam2hints --intronsonly`` on a coordinate-sorted BAM.

    Returns ``{chrom: {(start, end): mult}}``; the dict keeps the order in
    which the chromosomes first occur in the BAM, which is bam2hints'
    output order. Every record counts, as in bam2hints (secondary and
    supplementary alignments, both mates).
    """
    import pysam  # extras "align"

    out: Dict[str, Dict[Tuple[int, int], int]] = {}
    with pysam.AlignmentFile(str(bam_path), "rb", threads=max(1, threads)) as bam:
        names = bam.references
        for read in bam.fetch(until_eof=True):
            cig = read.cigartuples
            if not cig or read.reference_id < 0:
                continue
            # an intron needs a reference gap longer than MAX_GAP_LEN
            if not any(op in _REF_GAP and ln > MAX_GAP_LEN for op, ln in cig):
                continue
            introns = _alignment_introns(read.reference_start, cig)
            if not introns:
                continue
            d = out.setdefault(names[read.reference_id], {})
            for k in introns:
                d[k] = d.get(k, 0) + 1
    return out


def write_hints(bam_path: Path, genome: Path, out_path: Path, threads: int = 1) -> int:
    """Write ``hints.gff``; return the number of hint lines.

    Equivalent to ``bam2hints --intronsonly --in=BAM --out=tmp`` and
    ``filterIntronsFindStrand.pl genome tmp --score``: strand from the
    splice-site dinucleotides (GT-AG, GC-AG, AT-AC), other introns dropped,
    score column = multiplicity.
    """
    from varus.strand import StrandAssigner

    introns = bam2hints_introns(bam_path, threads=threads)
    strander = StrandAssigner(genome)
    tmp = out_path.with_name(out_path.name + ".tmp")
    n = n_drop = 0
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            for chrom, d in introns.items():
                for (s, e) in sorted(d):
                    strand = strander.resolve(chrom, s, e)
                    if strand is None:
                        n_drop += 1
                        continue
                    mult = d[(s, e)]
                    attr = (f"mult={mult};" if mult > 1 else "") + f"pri={PRIORITY};src={SOURCE}"
                    f.write(f"{chrom}\tb2h\tintron\t{s}\t{e}\t{mult}\t{strand}\t.\t{attr}\n")
                    n += 1
    finally:
        strander.close()
    os.replace(tmp, out_path)
    log.info("Intron hints: %s (%d introns; %d without a canonical splice site dropped)",
             out_path, n, n_drop)
    return n


# ---------------------------------------------------------------------------
# StringTie
# ---------------------------------------------------------------------------

def stringtie_cmd(out_gtf: Path, *, short_bam: Optional[Path] = None,
                  long_bam: Optional[Path] = None, threads: int = 1,
                  stringtie: str = "stringtie") -> List[str]:
    """The StringTie command line Paludamentum uses for this input."""
    cmd = [stringtie, "-p", str(max(1, threads)), "-o", str(out_gtf)]
    if short_bam and long_bam:
        return cmd + ["--mix", str(short_bam), str(long_bam)]   # long reads second
    if long_bam:
        return cmd + ["-L", str(long_bam)]
    if short_bam:
        return cmd + [str(short_bam)]
    raise ValueError("stringtie_cmd: no BAM given")


def run_stringtie(out_gtf: Path, *, short_bam: Optional[Path] = None,
                  long_bam: Optional[Path] = None, threads: int = 1,
                  stringtie: str = "stringtie") -> subprocess.Popen:
    """Start StringTie into ``out_gtf.tmp``; :func:`finish_stringtie` waits for it."""
    tmp = out_gtf.with_name(out_gtf.name + ".tmp")
    cmd = stringtie_cmd(tmp, short_bam=short_bam, long_bam=long_bam,
                        threads=threads, stringtie=stringtie)
    log.info("StringTie: %s", " ".join(cmd))
    return subprocess.Popen(cmd)


def finish_stringtie(proc: subprocess.Popen, out_gtf: Path) -> None:
    tmp = out_gtf.with_name(out_gtf.name + ".tmp")
    rc = proc.wait()
    if rc != 0 or not tmp.is_file():
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"stringtie exited with status {rc}")
    os.replace(tmp, out_gtf)
    n = sum(1 for line in open(out_gtf, encoding="utf-8")
            if "\ttranscript\t" in line)
    if n == 0:
        log.warning("StringTie assembled no transcript: %s", out_gtf)
    else:
        log.info("StringTie assembly: %s (%d transcripts)", out_gtf, n)


def content_md5(path: Path) -> str:
    """MD5 of the lines not starting with ``#`` (``hints.gff``)."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
    return h.hexdigest()


_GTF_VALUES = re.compile(r'(cov|longcov|FPKM|TPM) "([^"]+)"')
_GTF_TID = re.compile(r'transcript_id "([^"]+)"')


def gtf_fingerprint(path: Path) -> str:
    """MD5 of the assembled transcripts, independent of StringTie's IDs and order.

    StringTie with ``-p`` > 1 numbers genes and orders the output by thread
    completion, so the same BAM gives a different file each time (measured on
    brain, 2026-10-03: same transcripts and coverage values, other
    ``STRG`` ids and order). The fingerprint is over the sorted transcripts:
    sequence, strand, exons and the ``cov``/``longcov``/``FPKM``/``TPM`` values.
    """
    tx: Dict[str, list] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("transcript", "exon"):
                continue
            m = _GTF_TID.search(c[8])
            if not m:
                continue
            t = tx.setdefault(m.group(1), [c[0], c[6], [], ""])
            if c[2] == "exon":
                t[2].append(f"{c[3]}-{c[4]}")
            else:
                t[3] = ";".join(f"{k}={v}" for k, v in sorted(_GTF_VALUES.findall(c[8])))
    rows = sorted(f"{t[0]}\t{t[1]}\t{','.join(sorted(t[2], key=lambda e: int(e.split('-')[0])))}\t{t[3]}"
                  for t in tx.values())
    return hashlib.md5("\n".join(rows).encode()).hexdigest()


def assemble_bam(bam: Path, genome: Path, outdir: Path, *, longreads: bool,
                 threads: int = 1, stringtie: str = "stringtie") -> Dict[str, str]:
    """Write ``stringtie.gtf`` and ``hints.gff`` for one VARUS BAM.

    StringTie runs in the background while the hints are extracted. Returns
    the manifest keys that describe both files. Raises ``RuntimeError`` when
    StringTie fails; the BAM is left alone.
    """
    gtf = outdir / STRINGTIE_NAME
    hints = outdir / HINTS_NAME
    kw = {"long_bam": bam} if longreads else {"short_bam": bam}
    proc = run_stringtie(gtf, threads=threads, stringtie=stringtie, **kw)
    try:
        write_hints(bam, genome, hints, threads=max(1, min(4, threads)))
    except BaseException:
        proc.kill()
        proc.wait()
        gtf.with_name(gtf.name + ".tmp").unlink(missing_ok=True)
        raise
    finish_stringtie(proc, gtf)
    return assembly_header(gtf, hints, longreads=longreads, threads=threads,
                           stringtie=stringtie)


def assembly_header(gtf: Path, hints: Optional[Path], *, longreads: Optional[bool],
                    threads: int, stringtie: str = "stringtie") -> Dict[str, str]:
    """Manifest keys for an assembly; ``longreads=None`` means ``--mix``."""
    args = "--mix" if longreads is None else ("-L" if longreads else "")
    out = {
        "stringtie_version": tool_version(stringtie),
        "stringtie_args": args,
        "stringtie_threads": str(threads),
        "assembly": gtf.name,
        "assembly_md5": gtf_fingerprint(gtf),
    }
    if hints is not None:
        out["hints"] = hints.name
        out["hints_md5"] = content_md5(hints)
    return out


# ---------------------------------------------------------------------------
# `varus assemble`
# ---------------------------------------------------------------------------

ASSEMBLY_MANIFEST = "VARUS.assembly.tsv"


def _check_input(bam: Path, mode: str, genome_md5: str, skip_genome_check: bool) -> None:
    """The VARUS manifest next to ``bam`` must name this genome and mode."""
    from varus.provenance import MANIFEST_NAME, read_manifest

    man = bam.parent / MANIFEST_NAME
    if not man.is_file():
        log.warning("%s: no %s next to it; cannot check genome and mode", bam, MANIFEST_NAME)
        return
    header, _ = read_manifest(man)
    have = header.get("mode", "shortreads")
    if have != mode:
        raise SystemExit(f"{bam}: the VARUS run was in mode '{have}', expected '{mode}' "
                         f"(--short takes short-read runs, --long --longreads runs)")
    want = header.get("genome_md5", "")
    if want and not skip_genome_check and want != genome_md5:
        raise SystemExit(f"{bam} was aligned to another genome (MD5 {want}, here "
                         f"{genome_md5}). Pass the same FASTA; if it only differs in "
                         "formatting, pass --skip-genome-check.")


def assemble_cli(genome: Path, outdir: Path, *, short_bam: Optional[Path] = None,
                 long_bam: Optional[Path] = None, threads: int = 1,
                 skip_genome_check: bool = False, command: str = "") -> int:
    """``varus assemble``: one BAM → stringtie.gtf + hints.gff; both → ``--mix`` GTF."""
    from varus import __version__
    from varus.provenance import file_md5, now_iso

    if not short_bam and not long_bam:
        raise SystemExit("varus assemble: pass --short, --long or both")
    for b in (short_bam, long_bam):
        if b and not Path(b).is_file():
            raise SystemExit(f"{b}: not found")
    require_stringtie()
    md5 = file_md5(genome)
    if short_bam:
        _check_input(Path(short_bam), "shortreads", md5, skip_genome_check)
    if long_bam:
        _check_input(Path(long_bam), "longreads", md5, skip_genome_check)
    outdir.mkdir(parents=True, exist_ok=True)

    if short_bam and long_bam:
        gtf = outdir / STRINGTIE_NAME
        proc = run_stringtie(gtf, short_bam=short_bam, long_bam=long_bam, threads=threads)
        finish_stringtie(proc, gtf)
        keys = assembly_header(gtf, None, longreads=None, threads=threads)
        mode = "mixed"
    else:
        bam = Path(short_bam or long_bam)
        keys = assemble_bam(bam, genome, outdir, longreads=bool(long_bam), threads=threads)
        mode = "longreads" if long_bam else "shortreads"

    header = {
        "varus_assembly": "1",
        "varus_version": __version__,
        "created": now_iso(),
        "command": command,
        "genome": Path(genome).name,
        "genome_md5": md5,
        "mode": mode,
        "short_bam": str(short_bam or ""),
        "long_bam": str(long_bam or ""),
        **keys,
    }
    with open(outdir / ASSEMBLY_MANIFEST, "w", encoding="utf-8") as f:
        for k, v in header.items():
            v = str(v).replace("\n", " ").replace("\t", " ")
            f.write(f"#{k}={v}\n")
    log.info("Assembly written to %s (%s)", outdir, mode)
    return 0
