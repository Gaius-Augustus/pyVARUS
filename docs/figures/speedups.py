"""Draw docs/figures/speedups.svg: how the v2 speed-ups fit together.

    python docs/figures/speedups.py      # writes speedups.svg next to this file

Plain SVG, no dependencies. Numbers in the notes are from
docs/benchmark_logan.md (brain, 48 threads).
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape

W, H = 1200, 889

STYLE = {
    #          fill       stroke
    "main":  ("#dde8f8", "#3b6fb6"),
    "net":   ("#d7f0ea", "#2a8c74"),
    "align": ("#ebe2f8", "#6b4fa8"),
    "bg":    ("#ececec", "#777777"),
    "logan": ("#fde5c8", "#c9711c"),
    "file":  ("#ffffff", "#555555"),
}
NOTE = "#2e6b30"
out: list[str] = []


def box(x, y, w, h, lines, kind="main", size=12.5, bold_first=True, rx=6):
    fill, stroke = STYLE[kind]
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
               f'fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
    lh = size * 1.25
    y0 = y + h / 2 - (len(lines) - 1) * lh / 2 + size * 0.35
    for i, line in enumerate(lines):
        weight = ' font-weight="600"' if (i == 0 and bold_first) else ""
        out.append(f'<text x="{x + w / 2}" y="{y0 + i * lh:.1f}" font-size="{size}"'
                   f'{weight} text-anchor="middle">{escape(line)}</text>')


def text(x, y, s, size=12, anchor="start", color="#222", italic=False, weight=None):
    st = ' font-style="italic"' if italic else ""
    wt = f' font-weight="{weight}"' if weight else ""
    out.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" '
               f'text-anchor="{anchor}"{st}{wt}>{escape(s)}</text>')


def note(x, y, s, anchor="middle"):
    text(x, y, s, size=11, anchor=anchor, color=NOTE, italic=True)


def arrow(pts, color="#444", dash=False, width=1.5):
    d = "M " + " L ".join(f"{px},{py}" for px, py in pts)
    da = ' stroke-dasharray="5 4"' if dash else ""
    marker = "ahL" if color == STYLE["logan"][1] else "ah"
    out.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{da} '
               f'marker-end="url(#{marker})"/>')


def panel(y, h, title, subtitle):
    out.append(f'<rect x="12" y="{y}" width="{W - 24}" height="{h}" rx="10" '
               f'fill="none" stroke="#bbb" stroke-width="1.2"/>')
    text(28, y + 24, title, size=15, weight="700")
    text(28 + 9 * len(title) + 14, y + 24, subtitle, size=12, color="#555")


def lane(y, h, label, sub=""):
    out.append(f'<rect x="24" y="{y}" width="{W - 48}" height="{h}" rx="6" '
               f'fill="#fafafa" stroke="#e2e2e2"/>')
    text(34, y + 18, label, size=12, weight="600", color="#333")
    if sub:
        text(34, y + 33, sub, size=10.5, color="#666")


LOG = STYLE["logan"][1]

# ---------------------------------------------------------------- header
out.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
text(W / 2, 30, "pyVARUS v2: what runs concurrently, and where the time went",
     size=17, anchor="middle", weight="700")

# ---------------------------------------------------------------- panel A
panel(48, 212, "A   varus logan", "optional pre-screen, once per genome; 0.7–12 min in the benchmarks")
ya = 108
box(28, ya + 18, 96, 54, ["Runlist.tsv", "all SRA runs"], "file")
arrow([(124, ya + 45), (142, ya + 45)])
box(144, ya + 8, 128, 74, ["Sample ≤ 500", "runs by BioProject;", "HTTP HEAD (cached)"], "net")
arrow([(272, ya + 45), (290, ya + 45)])
box(292, ya + 8, 140, 74, ["Stream contigs", "8 S3 connections", "(~8 MB/s total)"], "net")
arrow([(432, ya + 45), (452, ya + 45)])
text(454, ya - 4, "chunks of 25 runs", size=10.5, color="#555")
for i in range(3):
    yy = ya + 4 + i * 30
    box(454, yy, 128, 25, [f"minimap2 group {i + 1}"], "align", size=11.5, bold_first=False, rx=4)
    arrow([(582, yy + 12.5), (606, yy + 12.5)])
    box(608, yy, 96, 25, ["scanner"], "main", size=11.5, bold_first=False, rx=4)
    arrow([(704, yy + 12.5), (714, yy + 12.5), (714, ya + 45), (728, ya + 45)] if i != 1
          else [(704, yy + 12.5), (728, yy + 12.5)], width=1.1)
note(579, ya + 108, "SAM piped, no BAM on disk; 3 groups: −20 %")
box(730, ya + 8, 146, 74, ["Gate", "tile breadth ≥ 10 %", "of the best run,", "divergence ≤ 0.05"], "main")
arrow([(876, ya + 45), (894, ya + 45)])
box(896, ya + 8, 124, 74, ["Greedy rank", "by the VARUS score", "S = Σ log(1 + cⱼ)"], "main")
arrow([(1020, ya + 45), (1036, ya + 45)])
for i, s in enumerate(["Runlist.logan.tsv", "per-run tile priors", "splice-DB seed"]):
    box(1038, ya - 2 + i * 32, 134, 27, [s], "logan", size=11.5, bold_first=False, rx=4)
note(875, ya + 108, "foreign runs dropped before any read download")
text(1105, ya + 108, "→ into B (orange)", size=11, anchor="middle", color=LOG, weight="600")
text(28, 242, "Threads: minimap2 gets --threads minus the scanners and one core for main + downloads; "
              "groups = 1 per ~15 minimap2 threads (1–4, fewer if the index copies do not fit in memory).",
     size=11, color="#444")

# ---------------------------------------------------------------- panel B
panel(274, 600, "B   varus run", "online loop, 1000 batches of 50 k reads; every stage below runs at the same time")

# lane 1: picks (main thread)
L1 = 316
lane(L1, 104, "Main thread", "choose")
box(160, L1 + 12, 200, 80, ["Estimator", "sparse per-run tile counts,", "λ = 3, a = 0.1"], "main")
box(302, L1 + 4, 62, 18, ["Logan prior"], "logan", size=9.5, bold_first=False, rx=9)
arrow([(360, L1 + 52), (384, L1 + 52)])
box(386, L1 + 12, 214, 80, ["Fresh pool", "untouched runs = 1 member;", "batch order built on first pick"], "main")
box(534, L1 + 4, 70, 18, ["Logan runlist"], "logan", size=9.5, bold_first=False, rx=9)
arrow([(600, L1 + 52), (624, L1 + 52)])
box(626, L1 + 12, 236, 80, ["Lazy greedy pick", "profit against observed counts", "+ expected gain of in-flight batches"], "main")
arrow([(862, L1 + 52), (886, L1 + 52)])
box(888, L1 + 12, 196, 80, ["Merge picks", "up to 10 consecutive batches", "of a run in one download"], "main")
note(260, L1 + 102, "Logan priors: 5.2 → 0.2 s per batch")
note(493, L1 + 102, "2 M mouse runs: 4 s → 0.02 s per pick")
note(940, L1 + 102, "loop −37 % (−63 % with Logan)", anchor="start")

# lane 2: downloads
L2 = 436
lane(L2, 84, "Download pool", "K = 6, one core")
for i in range(6):
    box(250 + i * 118, L2 + 18, 108, 44, ["fastq-dump", "spot range"], "net", size=11)
arrow([(930, L1 + 92), (930, L2 + 18)])
note(602, L2 + 78, "latency-bound (5–27 s per call); 6 in flight: 4.5–5× over serial")

# lane 3: aligner
L3 = 536
lane(L3, 84, "Aligner thread", "batch i + 1")
arrow([(360, L2 + 62), (360, L3 + 22)])
box(300, L3 + 16, 250, 48, ["HISAT2 --mm --no-unal", "-p = the budget below"], "align", size=11.5)
arrow([(550, L3 + 40), (574, L3 + 40)])
box(576, L3 + 16, 190, 48, ["samtools sort -@ ≤ 4", "BAM level 1"], "align", size=11.5)
note(980, L3 + 36, "aligns the next batch while the main thread")
note(980, L3 + 51, "counts this one (align-ahead): loop −23 %")

# lane 4: counting (main thread + scan workers)
L4 = 636
lane(L4, 118, "Main thread", "batch i")
arrow([(671, L3 + 64), (671, L3 + 90), (230, L3 + 90), (230, L4 + 22)])
box(160, L4 + 22, 150, 64, ["Quality gate", "≥ 5 % uniquely", "mapped reads"], "main", size=11.5)
arrow([(310, L4 + 54), (334, L4 + 54)])
box(336, L4 + 14, 250, 80, ["One-pass BAM scan", "UMRs, spliced reads, introns;", "merged batches: split by region", "over the scan workers"], "main", size=11.5)
arrow([(586, L4 + 54), (610, L4 + 54)])
box(612, L4 + 14, 230, 80, ["Splice-site DB", "strand cache; file rewritten", "only when new junctions appear"], "main", size=11.5)
box(776, L4 + 6, 70, 18, ["Logan seed"], "logan", size=9.5, bold_first=False, rx=9)
arrow([(842, L4 + 54), (866, L4 + 54)])
box(868, L4 + 22, 170, 64, ["Update tile counts", "and score S"], "main", size=11.5)
note(461, L4 + 110, "region split: scan 1.8–3× faster, same counts")
note(727, L4 + 110, "was 2 → 10 s per batch, growing")
# feedback: counts -> estimator; DB -> aligner
arrow([(1038, L4 + 54), (1160, L4 + 54), (1160, L1 - 8), (230, L1 - 8), (230, L1 + 12)], color="#3b6fb6")
text(1154, L3 + 10, "new counts", size=11, anchor="end", color="#3b6fb6", weight="600")
arrow([(727, L4 + 14), (727, L3 + 76), (425, L3 + 76), (425, L3 + 64)], dash=True, color="#6b4fa8")
text(733, L3 + 88, "DB for the next alignment", size=10.5, color="#6b4fa8")

# lane 5: merge
L5 = 770
lane(L5, 56, "Background", "merge")
box(160, L5 + 12, 280, 34, ["Rolling merge every 100 batches, -@ ≤ 4"], "bg", size=11.5, bold_first=False)
arrow([(440, L5 + 29), (464, L5 + 29)])
box(466, L5 + 12, 260, 34, ["Final merge of the parts, all threads"], "bg", size=11.5, bold_first=False)
arrow([(726, L5 + 29), (750, L5 + 29)])
box(752, L5 + 12, 100, 34, ["VARUS.bam"], "file", size=11.5)
arrow([(300, L4 + 86), (300, L5 + 12)], width=1.1)
text(306, L5 + 6, "accepted batch BAMs", size=10.5, color="#555")

# thread budget strip
yb = 836
text(34, yb + 16, "--threads 48:", size=11.5, weight="600")
x = 130
for n, label, kind in [(43, "HISAT2 -p 43", "align"), (4, "4 scan", "main"), (1, "", "net")]:
    w = n * 12.4
    box(x, yb + 3, w, 20, [label] if label else [], kind, size=10.5, bold_first=False, rx=3)
    x += w
text(x + 8, yb + 17, "+ 1 download core; a running rolling merge takes 4 from HISAT2", size=11, color="#444")

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">\n'
       '<defs>'
       '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
       'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#444"/></marker>'
       f'<marker id="ahL" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
       f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{LOG}"/></marker>'
       '</defs>\n' + "\n".join(out) + "\n</svg>\n")
Path(__file__).with_name("speedups.svg").write_text(svg, encoding="utf-8")
