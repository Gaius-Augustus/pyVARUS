"""CITATIONS.md: the references a run used, each with its DOI."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from tests.conftest import requires_pysam
from varus import citations
from varus.citations import CITATIONS_NAME, REFERENCES, citation_keys, write_citations


def test_every_reference_is_complete():
    for key, r in REFERENCES.items():
        assert re.fullmatch(r"10\.\d{4,9}/\S+", r["doi"]), key
        assert r["authors"] and r["title"] and r["journal"] and r["year"], key
        assert all(fam and giv for fam, giv in r["authors"]), key
        # a journal article has volume and pages; only a preprint may lack them
        assert (r["volume"] and r["pages"]) or r["journal"] == "bioRxiv", key
    dois = [r["doi"] for r in REFERENCES.values()]
    assert len(set(dois)) == len(dois)


@pytest.mark.parametrize("kw, expect", [
    (dict(mode="shortreads"),
     ["varus", "sra", "hisat2", "samtools", "stringtie"]),
    (dict(mode="longreads"),
     ["varus", "sra", "minimap2", "samtools", "stringtie"]),
    (dict(mode="mixed"),
     ["varus", "sra", "hisat2", "minimap2", "samtools", "stringtie"]),
    (dict(mode="shortreads", logan=True),
     ["varus", "sra", "hisat2", "minimap2", "samtools", "stringtie", "logan"]),
    (dict(mode="shortreads", assembled=False), ["varus", "sra", "hisat2", "samtools"]),
    (dict(mode="shortreads", extra=["nextflow", "nonsense"]),
     ["varus", "sra", "hisat2", "samtools", "stringtie", "nextflow"]),
])
def test_citation_keys(monkeypatch, kw, expect):
    monkeypatch.delenv(citations.EXTRA_ENV, raising=False)
    assert citation_keys(**kw) == expect


def test_citation_keys_from_environment(monkeypatch):
    monkeypatch.setenv(citations.EXTRA_ENV, "nextflow logan")
    keys = citation_keys(mode="shortreads")
    assert keys[-2:] == ["logan", "nextflow"] and "minimap2" in keys


def test_format_reference_and_bibtex():
    ref = citations.format_reference("stringtie")
    assert ref.startswith("- **StringTie:** Shinder I., Pertea G., Hu R., Rudnick Z., "
                          "Pertea M. (2026). ")
    assert ref.endswith("*Nature Methods*, 23(6):1126–1137. [DOI:10.1038/s41592-026-03080-3]"
                        "(https://doi.org/10.1038/s41592-026-03080-3)")
    assert citations._initials("Tsung-Cheng") == "T.-C."
    logan = citations.format_reference("logan")
    assert "Loll-Krippleber R., et al. (2024)" in logan and "*bioRxiv* (preprint)." in logan
    bib = citations.format_bibtex("logan")
    assert bib.startswith("@article{logan_2024,") and "Babaian, Artem" in bib
    assert "volume" not in bib and "doi = {10.1101/2024.07.30.605881}" in bib
    assert "pages = {907--915}" in citations.format_bibtex("hisat2")


def test_write_citations(tmp_path: Path, monkeypatch):
    monkeypatch.delenv(citations.EXTRA_ENV, raising=False)
    write_citations(tmp_path, mode="longreads", subcommand="replay")
    text = (tmp_path / CITATIONS_NAME).read_text(encoding="utf-8")
    assert "`varus replay`" in text and "```bibtex" in text
    keys = citation_keys(mode="longreads")
    for key in REFERENCES:
        doi = REFERENCES[key]["doi"]
        # twice in the list entry (text, link), twice in the BibTeX entry (doi, url)
        assert text.count(doi) == (4 if key in keys else 0), key
    # an unwritable directory is logged, not raised
    write_citations(tmp_path / "missing", mode="shortreads")


@requires_pysam
def test_run_writes_citations_and_logan_key(tmp_path: Path, monkeypatch):
    from tests.test_assemble import _run_controller

    monkeypatch.delenv(citations.EXTRA_ENV, raising=False)
    rc, out, header = _run_controller(tmp_path, monkeypatch)
    assert rc == 0 and header["logan"] == "0"
    text = (out / CITATIONS_NAME).read_text(encoding="utf-8")
    assert "10.1038/s41587-019-0201-4" in text and "10.1038/s41592-026-03080-3" in text
    assert "10.1101/2024.07.30.605881" not in text and "10.1093/bioinformatics/bty191" not in text
