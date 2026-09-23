"""Assign strand to intron hints using splice-site dinucleotides.

Reimplements ``filterIntronsFindStrand.pl`` (AUGUSTUS/BRAKER scripts).

For each intron at GFF coordinates (chrom, start, end) [1-based inclusive]:
- Donor dinucleotide: genome[start-1 : start+1]  (0-based, 2 chars)
- Acceptor dinucleotide: genome[end-2 : end]      (0-based, 2 chars)
- Concatenate → 4-char motif (lowercase)
- If motif ∈ allowed      → strand '+'
- If RC(motif) ∈ allowed  → strand '-'
- Otherwise               → intron is dropped (matching Perl behaviour)

The default allowed set matches the Perl default: gtag, gcag, atac.

Also provides :class:`StrandAssigner` (incremental, memoised variant used by
the controller) and :func:`write_hisat2_splice_sites` to convert an
:class:`~varus.introns.IntronCounts` into the tab-delimited format expected
by ``hisat2 --known-splicesite-infile``, and :func:`write_minimap2_junc_bed`
for the BED12 format expected by ``minimap2 --junc-bed``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import FrozenSet

from varus.introns import IntronCounts, IntronKey

log = logging.getLogger(__name__)

DEFAULT_ALLOWED: FrozenSet[str] = frozenset({"gtag", "gcag", "atac"})

_COMPLEMENT = str.maketrans("acgtACGT", "tgcaTGCA")


def _rc4(motif: str) -> str:
    """Reverse-complement a 4-char splice-site motif."""
    return motif[::-1].translate(_COMPLEMENT).lower()


class StrandAssigner:
    """Incremental, memoised splice-motif strand assignment.

    ``assign_strand`` re-slices the genome for *every* cumulative intron on
    each call, which made the per-batch intron-DB rebuild grow linearly with
    the number of introns seen (2 s -> 10 s per batch in production logs).
    This class opens the genome once and remembers the verdict for each
    ``(chrom, start, end)`` so the controller only pays for introns it has
    not seen before.
    """

    def __init__(
        self,
        genome_fasta: Path,
        allowed: FrozenSet[str] = DEFAULT_ALLOWED,
    ) -> None:
        self.genome_fasta = Path(genome_fasta)
        self.allowed = allowed
        self._fa = None
        # (chrom, start, end) -> "+" | "-" | None (None = non-canonical, dropped)
        self._cache: dict[tuple[str, int, int], str | None] = {}
        self._missing_chroms: set[str] = set()

    # -- genome access ---------------------------------------------------
    def _fasta(self):
        if self._fa is None:
            from pyfaidx import Fasta  # optional dependency (Linux/macOS)
            self._fa = Fasta(str(self.genome_fasta), as_raw=True)
        return self._fa

    def close(self) -> None:
        if self._fa is not None:
            try:
                self._fa.close()
            except Exception:  # pragma: no cover - best effort
                pass
            self._fa = None

    # -- cache management -------------------------------------------------
    @property
    def n_cached(self) -> int:
        return len(self._cache)

    def preload(self, stranded: IntronCounts) -> int:
        """Seed the cache from already-stranded introns (e.g. a Logan seed).

        Returns the number of keys added. Existing entries are kept.
        """
        n = 0
        for (chrom, start, end, strand), _ in stranded.counts.items():
            if strand not in ("+", "-"):
                continue
            key = (chrom, start, end)
            if key not in self._cache:
                self._cache[key] = strand
                n += 1
        return n

    def resolve(self, chrom: str, start: int, end: int) -> str | None:
        """Return "+", "-" or None (dropped) for a 1-based inclusive intron."""
        key = (chrom, start, end)
        try:
            return self._cache[key]
        except KeyError:
            pass
        fa = self._fasta()
        if chrom not in fa:
            if chrom not in self._missing_chroms:
                log.warning("Chrom %s absent from genome FASTA; dropping introns", chrom)
                self._missing_chroms.add(chrom)
            verdict: str | None = None
        else:
            donor = str(fa[chrom][start - 1 : start + 1]).lower()
            acceptor = str(fa[chrom][end - 2 : end]).lower()
            motif = donor + acceptor
            if motif in self.allowed:
                verdict = "+"
            elif _rc4(motif) in self.allowed:
                verdict = "-"
            else:
                verdict = None
        self._cache[key] = verdict
        return verdict

    def assign_new(self, introns: IntronCounts) -> tuple[IntronCounts, int]:
        """Strand a batch of introns; return (stranded, number of new keys).

        "New" means the ``(chrom, start, end)`` had not been resolved before
        this call. The caller uses the count to decide whether the aligner's
        splice-site DB (which carries no multiplicity) needs rewriting.
        """
        new: dict[IntronKey, int] = {}
        n_new = n_kept = n_dropped = 0
        for (chrom, start, end, _), mult in introns.counts.items():
            key = (chrom, start, end)
            seen = key in self._cache
            strand = self.resolve(chrom, start, end)
            if not seen:
                n_new += 1
            if strand is None:
                n_dropped += 1
                continue
            new[(chrom, start, end, strand)] = new.get((chrom, start, end, strand), 0) + mult
            n_kept += 1
        log.debug("assign_new: %d kept, %d dropped, %d new keys", n_kept, n_dropped, n_new)
        return IntronCounts(new), n_new


def assign_strand(
    introns: IntronCounts,
    genome_fasta: Path,
    allowed: FrozenSet[str] = DEFAULT_ALLOWED,
) -> IntronCounts:
    """Return a new IntronCounts with strand assigned; unrecognized introns dropped.

    Uses pyfaidx (optional dep) to fetch dinucleotides from the genome FASTA.
    Requires pyfaidx >= 0.7 and the FASTA to be indexed (``samtools faidx``
    or pyfaidx will build the index automatically on first run).

    Thin wrapper around :class:`StrandAssigner` for one-off use.
    """
    assigner = StrandAssigner(genome_fasta, allowed)
    try:
        stranded, _ = assigner.assign_new(introns)
    finally:
        assigner.close()
    n_kept = len(stranded.counts)
    n_dropped = len(introns.counts) - n_kept
    log.info("assign_strand: %d kept, %d dropped", n_kept, n_dropped)
    return stranded


def write_hisat2_splice_sites(introns: IntronCounts, path: Path) -> int:
    """Write a HISAT2 --known-splicesite-infile from stranded introns.

    Each line::

        chrom \\t donor_0based \\t acceptor_0based \\t strand

    where ``donor`` = last base of left exon (0-based) = start - 2,
    and ``acceptor`` = first base of right exon (0-based) = end.

    Introns with strand '.' are silently skipped.
    Returns the number of records written.
    """
    n = 0
    with path.open("w", encoding="utf-8") as f:
        for (chrom, start, end, strand), _ in sorted(introns.counts.items()):
            if strand not in ("+", "-"):
                continue
            donor = start - 2   # 0-based last base before the intron
            acceptor = end      # 0-based first base after the intron
            if donor < 0:
                continue
            f.write(f"{chrom}\t{donor}\t{acceptor}\t{strand}\n")
            n += 1
    log.info("Wrote %d splice sites to %s", n, path)
    return n


def write_minimap2_junc_bed(introns: IntronCounts, path: Path) -> int:
    """Write a BED12 junction file for ``minimap2 --junc-bed``.

    Each intron at 1-based inclusive coordinates ``(start, end)`` becomes one
    BED12 line representing two 1-bp anchors flanking the intron. minimap2
    reads the block boundaries to learn canonical splice sites.

    BED is 0-based half-open:

    - ``bed_start = start - 2`` (1 bp before the donor)
    - ``bed_end   = end + 1``   (1 bp after the acceptor)

    Two blocks of size 1 each, at offsets ``0`` and ``bed_end - bed_start - 1``.
    The intron multiplicity is used as the score (clipped to 1000).

    Introns with strand ``.`` or ``donor < 0`` are silently skipped (matches
    :func:`write_hisat2_splice_sites`). Returns the number of records written.
    """
    n = 0
    with path.open("w", encoding="utf-8") as f:
        for (chrom, start, end, strand), mult in sorted(introns.counts.items()):
            if strand not in ("+", "-"):
                continue
            bed_start = start - 2
            bed_end = end + 1
            if bed_start < 0:
                continue
            score = min(int(mult), 1000)
            block_count = 2
            block_sizes = "1,1"
            block_starts = f"0,{bed_end - bed_start - 1}"
            name = f"junc{n}"
            f.write(
                f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{strand}\t"
                f"{bed_start}\t{bed_end}\t0\t{block_count}\t{block_sizes}\t"
                f"{block_starts}\n"
            )
            n += 1
    log.info("Wrote %d junctions to %s", n, path)
    return n
