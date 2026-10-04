"""``CITATIONS.md``: what to cite for one pyVARUS run.

``varus run``, ``varus replay`` and ``varus assemble`` write it next to their
output. It lists only what that run used: HISAT2 or minimap2, StringTie,
Logan when the pre-screen was used, Nextflow when
the pipeline started the run.

Every entry of ``REFERENCES`` was checked against the Crossref record of its
DOI (2026-10-04). Add an entry only after checking it the same way.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, List

from varus import __version__

log = logging.getLogger(__name__)

CITATIONS_NAME = "CITATIONS.md"
REPO_URL = "https://github.com/Gaius-Augustus/pyVARUS"
# The Nextflow pipeline names here what `varus run` cannot see itself
# (space-separated keys of REFERENCES, e.g. "nextflow logan").
EXTRA_ENV = "VARUS_CITE_EXTRA"

# key -> label, authors (family, given), title, journal, year, volume, issue,
# pages, doi. `journal` "bioRxiv" marks a preprint (no volume/pages).
REFERENCES = {
    "varus": {
        "label": "VARUS",
        "authors": [("Stanke", "Mario"), ("Bruhn", "Willy"), ("Becker", "Felix"),
                    ("Hoff", "Katharina J.")],
        "title": "VARUS: sampling complementary RNA reads from the sequence read archive",
        "journal": "BMC Bioinformatics", "year": 2019, "volume": "20", "issue": "1",
        "pages": "558", "doi": "10.1186/s12859-019-3182-x",
    },
    "sra": {
        "label": "SRA",
        "authors": [("Katz", "Kenneth"), ("Shutov", "Oleg"), ("Lapoint", "Richard"),
                    ("Kimelman", "Michael"), ("Brister", "J. Rodney"),
                    ("O'Sullivan", "Christopher")],
        "title": "The Sequence Read Archive: a decade more of explosive growth",
        "journal": "Nucleic Acids Research", "year": 2022, "volume": "50", "issue": "D1",
        "pages": "D387-D390", "doi": "10.1093/nar/gkab1053",
    },
    "hisat2": {
        "label": "HISAT2",
        "authors": [("Kim", "Daehwan"), ("Paggi", "Joseph M."), ("Park", "Chanhee"),
                    ("Bennett", "Christopher"), ("Salzberg", "Steven L.")],
        "title": "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype",
        "journal": "Nature Biotechnology", "year": 2019, "volume": "37", "issue": "8",
        "pages": "907-915", "doi": "10.1038/s41587-019-0201-4",
    },
    "minimap2": {
        "label": "minimap2",
        "authors": [("Li", "Heng")],
        "title": "Minimap2: pairwise alignment for nucleotide sequences",
        "journal": "Bioinformatics", "year": 2018, "volume": "34", "issue": "18",
        "pages": "3094-3100", "doi": "10.1093/bioinformatics/bty191",
    },
    "samtools": {
        "label": "SAMtools",
        "authors": [("Danecek", "Petr"), ("Bonfield", "James K."), ("Liddle", "Jennifer"),
                    ("Marshall", "John"), ("Ohan", "Valeriu"), ("Pollard", "Martin O."),
                    ("Whitwham", "Andrew"), ("Keane", "Thomas"), ("McCarthy", "Shane A."),
                    ("Davies", "Robert M."), ("Li", "Heng")],
        "title": "Twelve years of SAMtools and BCFtools",
        "journal": "GigaScience", "year": 2021, "volume": "10", "issue": "2",
        "pages": "giab008", "doi": "10.1093/gigascience/giab008",
    },
    "stringtie": {
        "label": "StringTie",
        "authors": [("Shinder", "Ida"), ("Pertea", "Geo"), ("Hu", "Richard"),
                    ("Rudnick", "Zoe"), ("Pertea", "Mihaela")],
        "title": "StringTie3 improves total RNA-seq assembly by resolving nascent and "
                 "mature transcripts",
        "journal": "Nature Methods", "year": 2026, "volume": "23", "issue": "6",
        "pages": "1126-1137", "doi": "10.1038/s41592-026-03080-3",
    },
    "logan": {
        "label": "Logan",
        "authors": [("Chikhi", "Rayan"), ("Lemane", "Téo"), ("Loll-Krippleber", "Raphaël"),
                    ("Montoliu-Nerin", "Mercè"), ("Raffestin", "Brice"),
                    ("Camargo", "Antonio Pedro"), ("Miller", "Carson J."),
                    ("Fiamenghi", "Mateus Bernabe"), ("Agustinho", "Daniel Paiva"),
                    ("Majidian", "Sina"), ("Autric", "Greg"), ("Hugues", "Maxime"),
                    ("Lee", "Junkyoung"), ("Faure", "Roland"), ("Curry", "Kristen D."),
                    ("Moura de Sousa", "Jorge A."), ("Rocha", "Eduardo P. C."),
                    ("Koslicki", "David"), ("Medvedev", "Paul"), ("Gupta", "Purav"),
                    ("Shen", "Jessica"), ("Morales-Tapia", "Alejandro"), ("Sihuta", "Kate"),
                    ("Roy", "Peter J."), ("Brown", "Grant W."), ("Edgar", "Robert C."),
                    ("Korobeynikov", "Anton"), ("Steinegger", "Martin"),
                    ("Lareau", "Caleb A."), ("Peterlongo", "Pierre"), ("Babaian", "Artem")],
        "title": "Logan: Planetary-Scale Genome Assembly Surveys Life's Diversity",
        "journal": "bioRxiv", "year": 2024, "volume": "", "issue": "",
        "pages": "", "doi": "10.1101/2024.07.30.605881",
    },
    "nextflow": {
        "label": "Nextflow",
        "authors": [("Di Tommaso", "Paolo"), ("Chatzou", "Maria"), ("Floden", "Evan W."),
                    ("Prieto Barja", "Pablo"), ("Palumbo", "Emilio"),
                    ("Notredame", "Cedric")],
        "title": "Nextflow enables reproducible computational workflows",
        "journal": "Nature Biotechnology", "year": 2017, "volume": "35", "issue": "4",
        "pages": "316-319", "doi": "10.1038/nbt.3820",
    },
}

# Order of the entries in the file.
_ORDER = list(REFERENCES)
# More authors than this are cut to the first three and "et al." in the text
# list; the BibTeX entry always has all of them.
_MAX_AUTHORS = 11


def citation_keys(*, mode: str, logan: bool = False, assembled: bool = True,
                  extra: Iterable[str] = ()) -> List[str]:
    """Keys of ``REFERENCES`` for a run in ``mode`` (shortreads, longreads, mixed)."""
    keys = {"varus", "sra", "samtools"}
    if mode in ("shortreads", "mixed"):
        keys.add("hisat2")
    if mode in ("longreads", "mixed"):
        keys.add("minimap2")
    if assembled:
        keys.add("stringtie")
    if logan:
        # the pre-screen aligns Logan's contigs with minimap2
        keys.update(("logan", "minimap2"))
    keys.update(k for k in extra if k in REFERENCES)
    keys.update(k for k in os.environ.get(EXTRA_ENV, "").split() if k in REFERENCES)
    if "logan" in keys:
        keys.add("minimap2")
    return [k for k in _ORDER if k in keys]


def _initials(given: str) -> str:
    """'Joseph M.' -> 'J. M.', 'Tsung-Cheng' -> 'T.-C.'"""
    return " ".join("-".join(p[0] + "." for p in part.split("-") if p)
                    for part in given.split())


def format_reference(key: str) -> str:
    """One Markdown list entry: authors (year). title. journal, volume(issue):pages. DOI."""
    r = REFERENCES[key]
    authors = r["authors"]
    shown = authors if len(authors) <= _MAX_AUTHORS else authors[:3]
    names = ", ".join(f"{fam} {_initials(giv)}" for fam, giv in shown)
    if len(shown) < len(authors):
        names += ", et al."
    title = r["title"].rstrip(".")
    if r["volume"]:
        where = f"*{r['journal']}*, {r['volume']}({r['issue']}):{r['pages'].replace('-', '–')}."
    else:
        where = f"*{r['journal']}* (preprint)."
    return (f"- **{r['label']}:** {names} ({r['year']}). {title}. {where} "
            f"[DOI:{r['doi']}](https://doi.org/{r['doi']})")


def format_bibtex(key: str) -> str:
    r = REFERENCES[key]
    fields = [("author", " and ".join(f"{fam}, {giv}" for fam, giv in r["authors"])),
              ("title", "{" + r["title"] + "}"),
              ("journal", r["journal"]), ("year", str(r["year"])),
              ("volume", r["volume"]), ("number", r["issue"]),
              ("pages", r["pages"].replace("-", "--")), ("doi", r["doi"]),
              ("url", f"https://doi.org/{r['doi']}")]
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields if v)
    return f"@article{{{key}_{r['year']},\n{body}\n}}"


def write_citations(outdir: Path, *, mode: str, logan: bool = False,
                    assembled: bool = True, extra: Iterable[str] = (),
                    subcommand: str = "run") -> None:
    """Write ``CITATIONS.md`` to ``outdir``; a failure is logged, not raised."""
    keys = citation_keys(mode=mode, logan=logan, assembled=assembled, extra=extra)
    lines = [
        "# What to cite",
        "",
        f"Written by `varus {subcommand}` (pyVARUS {__version__}, {REPO_URL}). The list "
        "holds the method, the data sources and the tools that this run used. Please "
        "cite all of them, and name the pyVARUS version.",
        "",
        *(format_reference(k) for k in keys),
        "",
        "## BibTeX",
        "",
        "```bibtex",
        "\n\n".join(format_bibtex(k) for k in keys),
        "```",
        "",
    ]
    path = Path(outdir) / CITATIONS_NAME
    try:
        path.write_text("\n".join(lines), encoding="utf-8")
        log.info("What to cite for this run: %s", path)
    except OSError as e:
        log.error("Could not write %s: %s", path, e)
