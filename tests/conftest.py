"""Shared test fixtures and skip markers.

``requires_pysam`` skips a test if pysam is not importable. pysam needs htslib,
which doesn't build on stock Windows. Tests that build a BAM are gated on this
so the local Windows test pass stays green and the BAM tests still run on the
HPC/Linux CI.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


requires_pysam = pytest.mark.skipif(
    not _has("pysam"),
    reason="pysam not installed (extras 'align' on Linux/macOS)",
)

requires_zstandard = pytest.mark.skipif(
    not _has("zstandard"),
    reason="zstandard not installed (extras 'logan')",
)


def _fake_assemble_bam(bam, genome, outdir, *, longreads, threads=1, **kw):
    """Stand-in for varus.assemble.assemble_bam: no StringTie in most tests."""
    outdir = Path(outdir)
    (outdir / "stringtie.gtf").write_text("# fake\n")
    (outdir / "hints.gff").write_text("")
    return {"stringtie_version": "fake", "stringtie_args": "-L" if longreads else "",
            "stringtie_threads": str(threads), "assembly": "stringtie.gtf",
            "assembly_md5": "0", "hints": "hints.gff", "hints_md5": "0"}


@pytest.fixture(autouse=True)
def _no_stringtie(request, monkeypatch):
    """`varus run` / `replay` assemble VARUS.bam with StringTie at the end;
    the tests mock that unless marked ``real_assembly``."""
    if request.node.get_closest_marker("real_assembly"):
        return
    monkeypatch.setattr("varus.assemble.require_stringtie", lambda *a, **k: None)
    monkeypatch.setattr("varus.controller.assemble_bam", _fake_assemble_bam)
    monkeypatch.setattr("varus.replay.assemble_bam", _fake_assemble_bam)
