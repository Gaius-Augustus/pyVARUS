"""Tests for varus.assemble: bam2hints port, hints.gff, StringTie, `varus assemble`.

The comparisons with the original tools (``bam2hints``,
``filterIntronsFindStrand.pl``, ``stringtie``) run when they are found:
on PATH, or the Perl script via ``$FILTER_INTRONS_FIND_STRAND``.
"""

from __future__ import annotations

import os
import random
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from varus import assemble
from varus.assemble import (
    _alignment_introns, assemble_cli, bam2hints_introns, content_md5, gtf_fingerprint,
    stringtie_cmd, write_hints,
)
from varus.provenance import MANIFEST_NAME, file_md5, read_manifest, write_manifest

from tests.conftest import requires_pysam

real_assembly = pytest.mark.real_assembly

BAM2HINTS = shutil.which("bam2hints")
FILTER_PL = os.environ.get("FILTER_INTRONS_FIND_STRAND") or shutil.which(
    "filterIntronsFindStrand.pl")
STRINGTIE = shutil.which("stringtie")


def _c(s: str):
    """'10M5N10M' -> pysam cigartuples."""
    ops = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    out, num = [], ""
    for ch in s:
        if ch.isdigit():
            num += ch
        else:
            out.append((ops[ch], int(num)))
            num = ""
    return out


# ---------------------------------------------------------------------------
# bam2hints rules
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cigar,expect", [
    ("20M100N20M", [(121, 220)]),                    # start 100 (0-based)
    ("7M100N20M", []),                               # first block < 8
    ("20M100N7M", []),                               # last block < 8
    ("8M100N8M", [(109, 208)]),
    ("20M100N3M100N20M", [(121, 220), (224, 323)]),  # short inner block is fine
    ("20M10D20M100N20M", [(151, 250)]),              # deletion <= 14 closed
    ("20M14N20M", []),                               # short N is a deletion too
    ("20M20D20M100N20M", []),                        # gap 15..31 drops the alignment
    ("20M31N20M", []),
    ("20M32N20M", [(121, 152)]),                     # >= 32 is an intron, D or N
    ("20M40D20M", [(121, 160)]),
    ("20M350000N20M", [(121, 350120)]),
    ("20M350001N20M", []),                           # longer than maxintronlen
    ("20M6I100N20M", []),                            # query gap > 5 before the intron
    ("20M5I100N20M", [(121, 220)]),
    ("5S20M100N20M3H", [(121, 220)]),
    ("10=10X100N20M", [(121, 220)]),                 # =/X join to one block
    ("10M2I10M100N20M", [(121, 220)]),               # I splits a block, gap 0 joins it
    ("20M200000N20M200000N20M", []),                 # spans > maxgenelen
    ("40M", []),
])
def test_alignment_introns(cigar, expect):
    assert _alignment_introns(100, _c(cigar)) == expect


def _bam(path: Path, chroms, reads):
    """reads: (chrom_index, start0, cigar, flag[, tags]); written sorted."""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": ln} for n, ln in chroms]}
    reads = sorted(reads, key=lambda r: (r[0], r[1]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, r in enumerate(reads):
            ci, start, cig, flag = r[:4]
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.reference_id = ci
            a.reference_start = start
            a.cigarstring = cig
            qlen = sum(ln for op, ln in _c(cig) if op in (0, 1, 4, 7, 8))
            # not poly-A: StringTie -L trims poly-A ends of long reads
            a.query_sequence = "".join(random.Random(i).choices("ACGT", k=qlen))
            a.flag = flag
            a.mapping_quality = 60
            if len(r) > 4:
                a.set_tags(r[4])
            out.write(a)
    return path


def _random_cigar(rng: random.Random) -> str:
    parts = []
    if rng.random() < 0.2:
        parts.append(f"{rng.randint(1, 5)}{rng.choice('SH')}")
    for b in range(rng.randint(1, 5)):
        if b:
            kind = rng.random()
            if kind < 0.45:
                parts.append(f"{rng.choice([rng.randint(1, 14), rng.randint(15, 31), rng.randint(32, 400), rng.randint(32, 400), 350001])}N")
            elif kind < 0.6:
                parts.append(f"{rng.choice([rng.randint(1, 14), rng.randint(15, 60)])}D")
            elif kind < 0.75:
                parts.append(f"{rng.randint(1, 8)}I")
        n = rng.randint(1, 30)
        if rng.random() < 0.1:
            parts.append(f"{n}=")
            parts.append(f"{rng.randint(1, 3)}X")
        else:
            parts.append(f"{n}M")
    if rng.random() < 0.2:
        parts.append(f"{rng.randint(1, 5)}S")
    return "".join(parts)


@requires_pysam
@pytest.mark.skipif(BAM2HINTS is None, reason="bam2hints (AUGUSTUS) not on PATH")
def test_bam2hints_port_matches_bam2hints(tmp_path: Path):
    """Random CIGARs: same intron hints and multiplicities as the C++ tool."""
    rng = random.Random(7)
    templates = [_random_cigar(rng) for _ in range(400)]
    starts = [rng.randrange(0, 5000) for _ in range(25)]
    reads = [(rng.randrange(2), rng.choice(starts), rng.choice(templates),
              rng.choice([0, 16, 256, 2048])) for _ in range(4000)]
    bam = _bam(tmp_path / "a.bam", [("c1", 1_000_000), ("c2", 1_000_000)], reads)
    subprocess.run([BAM2HINTS, "--intronsonly", f"--in={bam}",
                    f"--out={tmp_path / 'b2h.gff'}"], check=True, capture_output=True)
    want = (tmp_path / "b2h.gff").read_text().splitlines()

    got = []
    for chrom, d in bam2hints_introns(bam).items():
        for (s, e) in sorted(d):
            m = d[(s, e)]
            got.append(f"{chrom}\tb2h\tintron\t{s}\t{e}\t0\t.\t.\t"
                       + (f"mult={m};" if m > 1 else "") + "pri=4;src=E")
    assert len(want) > 50
    assert got == want


# ---------------------------------------------------------------------------
# hints.gff (with strand)
# ---------------------------------------------------------------------------

def _genome(path: Path, seqs: dict) -> Path:
    path.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    return path


@requires_pysam
def test_write_hints_strand_and_score(tmp_path: Path):
    seq = list("c" * 600)
    seq[120:122] = "gt"; seq[218:220] = "ag"     # 121..220 GT-AG -> +
    seq[320:322] = "ct"; seq[418:420] = "ac"     # 321..420 CT-AC -> -
    genome = _genome(tmp_path / "g.fa", {"c1": "".join(seq)})   # 521..620 stays c..c
    reads = ([(0, 100, "20M100N20M", 0)] * 3 + [(0, 300, "20M100N20M", 0)]
             + [(0, 500, "20M100N20M", 0)] * 2)
    bam = _bam(tmp_path / "a.bam", [("c1", 700)], reads)
    n = write_hints(bam, genome, tmp_path / "hints.gff")
    assert n == 2
    assert (tmp_path / "hints.gff").read_text().splitlines() == [
        "c1\tb2h\tintron\t121\t220\t3\t+\t.\tmult=3;pri=4;src=E",
        "c1\tb2h\tintron\t321\t420\t1\t-\t.\tpri=4;src=E",
    ]


@requires_pysam
@pytest.mark.skipif(BAM2HINTS is None or FILTER_PL is None,
                    reason="bam2hints or filterIntronsFindStrand.pl not available")
def test_write_hints_matches_bam2hints_and_filter(tmp_path: Path):
    """hints.gff equals bam2hints --intronsonly | filterIntronsFindStrand.pl --score."""
    rng = random.Random(3)
    seqs = {n: [rng.choice("acgtACGT") for _ in range(ln)] for n, ln in (("c1", 8000), ("c2", 6000))}
    reads = []
    motifs = ["gt", "ag"], ["gc", "ag"], ["at", "ac"], ["ct", "ac"], ["ct", "gc"], ["gt", "at"]
    for ci, name in enumerate(seqs):
        for _ in range(60):                       # planted canonical introns
            s = rng.randrange(100, len(seqs[name]) - 700)
            ln = rng.randrange(40, 500)
            d, a = rng.choice(motifs)
            seqs[name][s:s + 2] = d.upper() if rng.random() < 0.5 else d
            seqs[name][s + ln - 2:s + ln] = a
            for _ in range(rng.randint(1, 4)):
                reads.append((ci, s - 30, f"30M{ln}N25M", 0))
        for _ in range(300):                      # random ones, mostly non-canonical
            reads.append((ci, rng.randrange(0, len(seqs[name]) - 1500),
                          _random_cigar(rng), 0))
    genome = _genome(tmp_path / "g.fa", {n: "".join(s) for n, s in seqs.items()})
    bam = _bam(tmp_path / "a.bam", [(n, len(s)) for n, s in seqs.items()], reads)
    subprocess.run([BAM2HINTS, "--intronsonly", f"--in={bam}",
                    f"--out={tmp_path / 'b2h.gff'}"], check=True, capture_output=True)
    want = subprocess.run(["perl", FILTER_PL, str(genome), str(tmp_path / "b2h.gff"), "--score"],
                          check=True, capture_output=True, text=True).stdout.splitlines()
    write_hints(bam, genome, tmp_path / "hints.gff")
    got = (tmp_path / "hints.gff").read_text().splitlines()
    assert len(want) > 60
    assert got == want


# ---------------------------------------------------------------------------
# StringTie
# ---------------------------------------------------------------------------

def test_stringtie_cmd():
    o = Path("o.gtf")
    assert stringtie_cmd(o, short_bam=Path("s.bam"), threads=8) == [
        "stringtie", "-p", "8", "-o", "o.gtf", "s.bam"]
    assert stringtie_cmd(o, long_bam=Path("l.bam")) == [
        "stringtie", "-p", "1", "-o", "o.gtf", "-L", "l.bam"]
    assert stringtie_cmd(o, short_bam=Path("s.bam"), long_bam=Path("l.bam"), threads=2) == [
        "stringtie", "-p", "2", "-o", "o.gtf", "--mix", "s.bam", "l.bam"]
    with pytest.raises(ValueError):
        stringtie_cmd(o)


def _fake_stringtie(tmp_path: Path, body: str) -> str:
    exe = tmp_path / "fake_stringtie"
    exe.write_text("#!/bin/sh\n" + body)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return str(exe)


_FAKE_OK = """out=""
while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out="$2"; shift; done
printf '# stringtie %s\\nc1\\tStringTie\\ttranscript\\t1\\t9\\t1000\\t+\\t.\\tgene_id "STRG.1";\\n' "$out" > "$out"
"""


def test_finish_stringtie_renames_and_fails_cleanly(tmp_path: Path):
    gtf = tmp_path / "stringtie.gtf"
    p = assemble.run_stringtie(gtf, short_bam=tmp_path / "x.bam",
                               stringtie=_fake_stringtie(tmp_path, _FAKE_OK))
    assemble.finish_stringtie(p, gtf)
    assert gtf.is_file() and not gtf.with_name("stringtie.gtf.tmp").exists()

    gtf2 = tmp_path / "b.gtf"
    p = assemble.run_stringtie(gtf2, short_bam=tmp_path / "x.bam",
                               stringtie=_fake_stringtie(tmp_path, _FAKE_OK + "exit 1\n"))
    with pytest.raises(RuntimeError, match="status 1"):
        assemble.finish_stringtie(p, gtf2)
    assert not gtf2.exists() and not gtf2.with_name("b.gtf.tmp").exists()


def test_gtf_fingerprint_ignores_ids_and_order(tmp_path: Path):
    def tx(tid, s, e, cov, exons):
        out = [f'c1\tStringTie\ttranscript\t{s}\t{e}\t1000\t+\t.\tgene_id "G"; '
               f'transcript_id "{tid}"; cov "{cov}"; FPKM "1.0"; TPM "2.0";']
        out += [f'c1\tStringTie\texon\t{a}\t{b}\t1000\t+\t.\tgene_id "G"; '
                f'transcript_id "{tid}"; exon_number "1"; cov "{cov}";' for a, b in exons]
        return out
    t1 = tx("STRG.1.1", 10, 90, "5.0", [(10, 20), (50, 90)])
    t2 = tx("STRG.2.1", 200, 300, "7.0", [(200, 300)])
    a, b, c = tmp_path / "a", tmp_path / "b", tmp_path / "c"
    a.write_text("# stringtie -o /x/a\n" + "\n".join(t1 + t2) + "\n")
    b.write_text("# stringtie -o /y/b\n" + "\n".join(
        [l.replace("STRG.2.1", "STRG.9.1") for l in t2] +
        [l.replace("STRG.1.1", "STRG.4.1") for l in t1[:1] + t1[:0:-1]]) + "\n")
    assert gtf_fingerprint(a) == gtf_fingerprint(b)
    c.write_text("\n".join(t1 + [l.replace('"7.0"', '"7.5"') for l in t2]) + "\n")
    assert gtf_fingerprint(a) != gtf_fingerprint(c)


def test_content_md5_ignores_comments(tmp_path: Path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("# stringtie -o /x/a\nline\n")
    b.write_text("# stringtie -o /y/b\nline\n")
    assert content_md5(a) == content_md5(b)
    b.write_text("# x\nother\n")
    assert content_md5(a) != content_md5(b)


@requires_pysam
@real_assembly
def test_assemble_bam_writes_both(tmp_path: Path, monkeypatch):
    seq = list("c" * 600)
    seq[120:122] = "gt"; seq[218:220] = "ag"
    genome = _genome(tmp_path / "g.fa", {"c1": "".join(seq)})
    bam = _bam(tmp_path / "VARUS.bam", [("c1", 700)], [(0, 100, "20M100N20M", 0)])
    keys = assemble.assemble_bam(bam, genome, tmp_path, longreads=True, threads=3,
                                 stringtie=_fake_stringtie(tmp_path, _FAKE_OK))
    assert (tmp_path / "stringtie.gtf").is_file() and (tmp_path / "hints.gff").is_file()
    assert keys["stringtie_args"] == "-L" and keys["stringtie_threads"] == "3"
    assert keys["assembly_md5"] == gtf_fingerprint(tmp_path / "stringtie.gtf")
    assert keys["hints_md5"] == content_md5(tmp_path / "hints.gff")


def _transcript_reads(rng, exons, n, length, strand="+"):
    """Reads of ``length`` drawn along a spliced transcript (0-based half-open exons)."""
    tx = [(s, e) for s, e in exons]
    tlen = sum(e - s for s, e in tx)
    out = []
    for _ in range(n):
        off = rng.randrange(0, tlen - length + 1)
        pos, cig, left, start = 0, [], length, None
        for i, (s, e) in enumerate(tx):
            if off >= pos + (e - s):
                pos += e - s
                continue
            a = s + max(0, off - pos)
            if start is None:
                start = a
            elif cig:
                cig.append(f"{a - prev_end}N")
            take = min(e - a, left)
            cig.append(f"{take}M")
            left -= take
            prev_end = a + take
            pos += e - s
            if left == 0:
                break
        out.append((0, start, "".join(cig), 0, [("XS", strand)]))
    return out


def _gtf_introns(path: Path):
    """Intron chain of the assembled transcripts (StringTie trims low-coverage ends)."""
    exons = sorted((int(f[3]), int(f[4])) for f in (
        line.split("\t") for line in path.read_text().splitlines() if not line.startswith("#"))
        if f[2] == "exon")
    return [(e1 + 1, s2 - 1) for (_, e1), (s2, _) in zip(exons, exons[1:])]


@requires_pysam
@real_assembly
@pytest.mark.skipif(STRINGTIE is None, reason="stringtie not on PATH")
def test_real_stringtie_short_long_mix(tmp_path: Path):
    rng = random.Random(1)
    seq = [rng.choice("acgt") for _ in range(5000)]
    exons = [(1000, 1200), (1500, 1700), (2000, 2300)]
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        seq[e1:e1 + 2] = "gt"
        seq[s2 - 2:s2] = "ag"
    genome = _genome(tmp_path / "g.fa", {"c1": "".join(seq)})
    short = _bam(tmp_path / "short.bam", [("c1", 5000)], _transcript_reads(rng, exons, 400, 100))
    long = _bam(tmp_path / "long.bam", [("c1", 5000)], _transcript_reads(rng, exons, 60, 500))
    want = [(e1 + 1, s2) for (_, e1), (s2, _) in zip(exons, exons[1:])]

    for name, kw in (("s", {"longreads": False}), ("l", {"longreads": True})):
        d = tmp_path / name
        d.mkdir()
        bam = long if kw["longreads"] else short
        keys = assemble.assemble_bam(bam, genome, d, threads=2, **kw)
        assert _gtf_introns(d / "stringtie.gtf") == want, name
        assert keys["stringtie_version"] != "not found"
        assert (d / "hints.gff").read_text().count("\t+\t") == 2

    d = tmp_path / "m"
    d.mkdir()
    p = assemble.run_stringtie(d / "stringtie.gtf", short_bam=short, long_bam=long)
    assemble.finish_stringtie(p, d / "stringtie.gtf")
    assert _gtf_introns(d / "stringtie.gtf") == want


# ---------------------------------------------------------------------------
# `varus assemble`
# ---------------------------------------------------------------------------

def _varus_dir(tmp_path: Path, name: str, genome: Path, mode: str) -> Path:
    d = tmp_path / name
    d.mkdir()
    (d / "VARUS.bam").write_bytes(b"BAM")
    write_manifest(d / MANIFEST_NAME, {"genome_md5": file_md5(genome), "mode": mode}, [])
    return d / "VARUS.bam"


def test_assemble_cli_checks_and_mix(tmp_path: Path, monkeypatch):
    genome = _genome(tmp_path / "g.fa", {"c1": "acgt" * 10})
    short = _varus_dir(tmp_path, "s", genome, "shortreads")
    long = _varus_dir(tmp_path, "l", genome, "longreads")
    started = []

    def fake_run(out_gtf, **kw):
        started.append(kw)
        out_gtf.write_text("# h\nc1\tStringTie\ttranscript\t1\t9\n")
        return None

    monkeypatch.setattr(assemble, "run_stringtie", fake_run)
    monkeypatch.setattr(assemble, "finish_stringtie", lambda p, g: None)

    with pytest.raises(SystemExit, match="mode 'longreads'"):
        assemble_cli(genome, tmp_path / "o", short_bam=long)
    other = _genome(tmp_path / "h.fa", {"c1": "aaaa"})
    with pytest.raises(SystemExit, match="another genome"):
        assemble_cli(other, tmp_path / "o", short_bam=short, long_bam=long)
    with pytest.raises(SystemExit, match="--short, --long"):
        assemble_cli(genome, tmp_path / "o")

    out = tmp_path / "mix"
    assert assemble_cli(genome, out, short_bam=short, long_bam=long, threads=3) == 0
    assert started == [{"short_bam": short, "long_bam": long, "threads": 3}]
    man = (out / assemble.ASSEMBLY_MANIFEST).read_text()
    assert f"#genome_md5={file_md5(genome)}" in man
    assert "#mode=mixed" in man and "#stringtie_args=--mix" in man
    assert "#hints=" not in man


# ---------------------------------------------------------------------------
# varus run: finalize with assembly / --drop-bam
# ---------------------------------------------------------------------------

def _run_controller(tmp_path: Path, monkeypatch, **kw):
    import random as _r
    from varus.controller import Controller, RunState
    from tests.test_speedups import _cfg, _rec, _write_genome
    from tests.test_provenance import _mock_growing_db

    tmp_path.mkdir(parents=True, exist_ok=True)
    cfg = _cfg(tmp_path, max_batches=3, **kw)
    _write_genome(cfg.genome)
    runs = [RunState.from_record(_rec(f"R{i}", spots=5_000_000), cfg.batch_size, _r.Random(0))
            for i in range(2)]
    _mock_growing_db(monkeypatch, n_runs=2)
    rc = Controller(cfg, runs).run()
    header, _ = read_manifest(cfg.outdir / MANIFEST_NAME)
    return rc, cfg.outdir, header


@requires_pysam
def test_run_assembles_and_drops_bam(tmp_path: Path, monkeypatch):
    rc, out, header = _run_controller(tmp_path / "keep", monkeypatch)
    assert rc == 0 and (out / "VARUS.bam").is_file() and (out / "stringtie.gtf").is_file()
    assert header["assembly"] == "stringtie.gtf" and header["hints"] == "hints.gff"
    assert "bam_dropped" not in header

    rc, out, header = _run_controller(tmp_path / "drop", monkeypatch, drop_bam=True)
    assert rc == 0 and not (out / "VARUS.bam").exists()
    assert (out / "stringtie.gtf").is_file() and header["bam_dropped"] == "1"


@requires_pysam
def test_run_assembly_failure_keeps_bam(tmp_path: Path, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("stringtie exited with status 1")

    monkeypatch.setattr("varus.controller.assemble_bam", boom)
    rc, out, header = _run_controller(tmp_path, monkeypatch, drop_bam=True)
    assert rc == assemble.EXIT_ASSEMBLY_FAILED
    assert (out / "VARUS.bam").is_file()
    assert "assembly" not in header and header["batches"]


@real_assembly
def test_require_stringtie(monkeypatch):
    monkeypatch.setattr(assemble.shutil, "which", lambda t: None)
    with pytest.raises(SystemExit, match="StringTie 3.0.3"):
        assemble.require_stringtie()


@real_assembly
def test_varus_run_fails_early_without_stringtie(tmp_path: Path, monkeypatch):
    from varus import cli

    (tmp_path / "g.fa").write_text(">c1\nacgt\n")
    monkeypatch.setattr(assemble.shutil, "which", lambda t: None)
    with pytest.raises(SystemExit, match="stringtie not found"):
        cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--runlist", str(tmp_path / "r"),
                  "--index", str(tmp_path / "i"), "--no-logan"])


def test_replay_assembles_with_manifest_threads_and_reports_difference(tmp_path: Path,
                                                                      monkeypatch, caplog):
    from varus import replay as rp
    from tests.conftest import _fake_assemble_bam

    seen = {}

    def rec(bam, genome, outdir, *, longreads, threads=1, **kw):
        seen.update(longreads=longreads, threads=threads)
        return _fake_assemble_bam(bam, genome, outdir, longreads=longreads, threads=threads)

    monkeypatch.setattr(rp, "assemble_bam", rec)
    cfg = rp.ReplayConfig(manifest=tmp_path / "m", genome=tmp_path / "g.fa",
                          index=tmp_path / "i", outdir=tmp_path, threads=2)
    header = {"stringtie_threads": "7", "stringtie_version": "fake",
              "assembly_md5": "1", "hints_md5": "0"}
    with caplog.at_level("WARNING"):
        assert rp._assemble(cfg, header, tmp_path / "VARUS.bam", True) == 0
    assert seen == {"longreads": True, "threads": 7}
    assert "stringtie.gtf differs" in caplog.text and "hints.gff differs" not in caplog.text

    monkeypatch.setattr(rp, "assemble_bam", lambda *a, **k: (_ for _ in ()).throw(
        RuntimeError("stringtie exited with status 1")))
    assert rp._assemble(cfg, {}, tmp_path / "VARUS.bam", False) == assemble.EXIT_ASSEMBLY_FAILED
