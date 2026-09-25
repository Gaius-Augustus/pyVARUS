"""Count uniquely-mapped reads per genome tile from a sorted BAM.

The tile is the proxy for a "transcribed unit" in the VARUS algorithm: the
genome is partitioned into non-overlapping windows of ``tile_size`` bp
(default 5 kb). For each read, all of its alignments are pooled by tile; the
read is counted iff:

* it mapped to exactly one tile, **and**
* that tile saw at most two alignment records for the read (so a paired-end
  pair both landing on the same tile still counts once, but a multi-mapper
  is excluded).

This matches ``RNAread::UMR()`` in the legacy code (see
``legacy/Implementation/src/RNAread.cpp``).

Coordinate convention
---------------------
The legacy reads the SAM POS column directly (1-based). pysam exposes
``reference_start`` as 0-based, so the equivalent tile index is
``(reference_start + 1) // tile_size``.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from varus.introns import IntronCounts, IntronKey, iter_introns

log = logging.getLogger(__name__)

Tile = Tuple[str, int]  # (chromosome, tile_index)


@dataclass
class BAMStats:
    """Summary statistics from a single BAM-pass."""
    umr_counts: Dict[Tile, int]
    n_reads: int           # distinct mapped read names
    n_spliced: int         # reads with at least one N (intron) op


def count_umrs_per_tile(bam_path: Path, tile_size: int) -> Dict[Tile, int]:
    """Return ``{(chrom, tile_idx): umr_count}`` for the given BAM.

    Reference API: the controller uses :func:`scan_batch_bam`; this and
    :func:`count_bam_stats` are kept for the equivalence tests and scripts.
    """
    return count_bam_stats(bam_path, tile_size).umr_counts


def count_bam_stats(bam_path: Path, tile_size: int) -> BAMStats:
    """Single-pass BAM scan: UMR counts per tile + spliced-read count.

    Reference API (see :func:`count_umrs_per_tile`).

    Returns a :class:`BAMStats` with:
    - ``umr_counts``: ``{(chrom, tile_idx): n}`` — uniquely-mapped read count.
    - ``n_reads``: total distinct mapped read names.
    - ``n_spliced``: reads with at least one CIGAR N-op (intron).
    """
    return scan_batch_bam(bam_path, tile_size)[0]


def scan_batch_bam(bam_path: Path, tile_size: int) -> Tuple[BAMStats, IntronCounts]:
    """One pass over a batch BAM: :class:`BAMStats` plus intron multiplicities.

    Same results as :func:`count_bam_stats` followed by
    :func:`varus.introns.extract_introns_from_bam`, which each read the whole
    BAM; the loop scanned every batch twice (about 0.5 s per 50 000 spots).
    Introns are counted from every record, secondary and supplementary
    included, as in ``extract_introns_from_bam``.
    """
    if tile_size <= 0:
        raise ValueError("tile_size must be > 0")

    import pysam  # local import: pysam is an optional install (extras "align")

    per_read: dict[str, dict[Tile, int]] = defaultdict(lambda: defaultdict(int))
    spliced_reads: set[str] = set()
    intron_counts: Dict[IntronKey, int] = {}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue
            spliced = False
            for chrom, start, end in iter_introns(read):
                ikey: IntronKey = (chrom, start, end, ".")
                intron_counts[ikey] = intron_counts.get(ikey, 0) + 1
                spliced = True
            tile = (read.reference_start + 1) // tile_size
            key: Tile = (read.reference_name, tile)
            per_read[read.query_name][key] += 1
            if spliced:
                spliced_reads.add(read.query_name)

    umr_counts: Dict[Tile, int] = defaultdict(int)
    for read_name, tiles in per_read.items():
        if len(tiles) != 1:
            continue
        ((tile, hits),) = tiles.items()
        if hits <= 2:  # one pair landing on the same tile is allowed
            umr_counts[tile] += 1

    n_reads = len(per_read)
    n_spliced = len(spliced_reads)
    log.info(
        "BAM %s: %d UMRs across %d tiles; %d/%d reads spliced",
        bam_path, sum(umr_counts.values()), len(umr_counts), n_spliced, n_reads,
    )
    log.info("BAM %s: %d distinct introns", bam_path, len(intron_counts))
    stats = BAMStats(
        umr_counts=dict(umr_counts),
        n_reads=n_reads,
        n_spliced=n_spliced,
    )
    return stats, IntronCounts(intron_counts)


# ---------------------------------------------------------------------------
# Parallel scan of large (merged) batch BAMs
# ---------------------------------------------------------------------------

Region = Tuple[str, int, int]  # (chrom, start, end), 0-based half-open


def split_regions(
    references: List[Tuple[str, int]], n_parts: int, tile_size: int
) -> List[List[Region]]:
    """Cut the genome into ``n_parts`` groups of regions of similar length.

    Region boundaries fall on tile multiples; small references are packed
    together. Used to distribute a coordinate-sorted, indexed BAM over
    worker processes.
    """
    total = sum(max(0, ln) for _, ln in references)
    n_parts = max(1, int(n_parts))
    if total == 0:
        return []
    target = max(tile_size, -(-total // n_parts))
    target = -(-target // tile_size) * tile_size
    groups: List[List[Region]] = [[]]
    fill = 0
    for chrom, ln in references:
        start = 0
        while start < ln:
            room = max(tile_size, (target - fill) // tile_size * tile_size)
            take = min(ln - start, room)
            groups[-1].append((chrom, start, start + take))
            fill += take
            start += take
            if fill >= target:
                groups.append([])
                fill = 0
    return [g for g in groups if g]


def _scan_regions(bam_path: str, regions: List[Region], tile_size: int):
    """Worker body: scan the records *starting* in ``regions``.

    Reads whose alignments are all here (HISAT2 ``NH:i:1`` and, for pairs
    with a mapped mate, both mates seen) are decided locally, exactly as in
    :func:`scan_batch_bam`. Multimappers and pairs split across regions are
    returned by name for the caller to merge.
    """
    import pysam

    per_read: Dict[str, list] = {}   # name -> [tiles, all NH==1, seen, expected, spliced]
    intron_counts: Dict[IntronKey, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, start, end in regions:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.reference_start < start:
                    continue
                spliced = False
                for c, s, e in iter_introns(read):
                    ikey: IntronKey = (c, s, e, ".")
                    intron_counts[ikey] = intron_counts.get(ikey, 0) + 1
                    spliced = True
                name = read.query_name
                rec = per_read.get(name)
                if rec is None:
                    expected = 2 if read.is_paired and not read.mate_is_unmapped else 1
                    rec = per_read[name] = [{}, True, 0, expected, False]
                key: Tile = (read.reference_name, (read.reference_start + 1) // tile_size)
                t = rec[0]
                t[key] = t.get(key, 0) + 1
                if not (read.has_tag("NH") and read.get_tag("NH") == 1):
                    rec[1] = False
                rec[2] += 1
                if spliced:
                    rec[4] = True

    umr: Dict[Tile, int] = {}
    n_reads = n_spliced = 0
    pending: Dict[str, Tuple[Dict[Tile, int], bool]] = {}
    for name, (t, nh1, seen, expected, spliced) in per_read.items():
        if nh1 and seen == expected:
            n_reads += 1
            n_spliced += spliced
            if len(t) == 1:
                ((tile, hits),) = t.items()
                if hits <= 2:
                    umr[tile] = umr.get(tile, 0) + 1
        else:
            pending[name] = (t, spliced)
    return umr, n_reads, n_spliced, pending, intron_counts


def scan_batch_bam_parallel(
    bam_path: Path, tile_size: int, executor, n_parts: int
) -> Tuple[BAMStats, IntronCounts]:
    """:func:`scan_batch_bam` split by genome region over ``executor``.

    Same result as :func:`scan_batch_bam` for HISAT2 BAMs (``NH`` tags); reads
    without an ``NH`` tag are simply merged centrally. The BAM must be
    coordinate-sorted; it is indexed here if no index exists.
    """
    if tile_size <= 0:
        raise ValueError("tile_size must be > 0")
    import pysam

    bam_path = Path(bam_path)
    if not (Path(str(bam_path) + ".bai").is_file() or Path(str(bam_path) + ".csi").is_file()):
        pysam.index(str(bam_path))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        refs = list(zip(bam.references, bam.lengths))
    groups = split_regions(refs, n_parts, tile_size)
    futures = [executor.submit(_scan_regions, str(bam_path), g, tile_size) for g in groups]

    umr: Dict[Tile, int] = defaultdict(int)
    intron_counts: Dict[IntronKey, int] = defaultdict(int)
    n_reads = n_spliced = 0
    merged: Dict[str, list] = {}
    for fut in futures:
        u, nr, ns, pending, ic = fut.result()
        for k, v in u.items():
            umr[k] += v
        for k, v in ic.items():
            intron_counts[k] += v
        n_reads += nr
        n_spliced += ns
        for name, (t, spliced) in pending.items():
            m = merged.get(name)
            if m is None:
                merged[name] = [dict(t), spliced]
            else:
                mt = m[0]
                for k, v in t.items():
                    mt[k] = mt.get(k, 0) + v
                m[1] = m[1] or spliced
    for t, spliced in merged.values():
        n_reads += 1
        n_spliced += spliced
        if len(t) == 1:
            ((tile, hits),) = t.items()
            if hits <= 2:
                umr[tile] += 1
    log.info(
        "BAM %s: %d UMRs across %d tiles; %d/%d reads spliced (%d workers, %d merged centrally)",
        bam_path, sum(umr.values()), len(umr), n_spliced, n_reads, len(groups), len(merged),
    )
    stats = BAMStats(umr_counts=dict(umr), n_reads=n_reads, n_spliced=n_spliced)
    return stats, IntronCounts(dict(intron_counts))
