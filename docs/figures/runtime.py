"""Draw docs/figures/runtime.svg: wall time of pyVARUS v2 on the benchmark genomes.

    python docs/figures/runtime.py      # writes runtime.svg next to this file

Numbers are from docs/benchmark_logan.md (brain, 1000 batches). Wall time is
`varus logan` plus `varus run` (including the final merge); bars are seed
means, whiskers the range over seeds.
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape

# (label, logan stage minutes, varus run minutes, (min total, max total) or None, note)
GENOMES = [
    ("Coelastrella tenuitheca", "9 runs", [
        ("v1 behaviour (A0)", 0, 272.4, None, ""),
        ("v2", 0, 16.1, (15.8, 16.6), ""),
        ("v2 + Logan", 0.7, 14.1, (14.7, 14.9), "same S (100 %)"),
    ]),
    ("Chlorella sorokiniana", "391 runs", [
        ("v1 behaviour (A0)", 0, 267.6, None, ""),
        ("v2", 0, 24.8, (23.1, 25.8), "S 99.6 %, 409 batches rejected"),
        ("v2 + Logan", 8.6, 20.3, (23.1, 32.1), "S 106 %, 12 rejected"),
    ]),
    ("Drosophila melanogaster", "115 k runs", [
        ("v1 behaviour (A0)", 0, 259.5, None, ""),
        ("v2", 0, 18.3, (16.5, 19.4), "S 77–89 %"),
        ("v2 + Logan", 11.5, 19.2, (29.7, 31.4), "S 101 %"),
    ]),
    ("Mus musculus", "2 M runs", [
        ("v2", 0, 17.3, (15.1, 19.9), "intron Sn 0.83–0.86"),
        ("v2 + Logan", 75.6, 18.6, (92.2, 96.0), "intron Sn 0.91"),
    ]),
]
# Chlorella sorokiniana, seed 1: (threads, logan stage, run with Logan, run without)
THREADS = [(4, 64.7, 37.1, 32.6), (8, 35.6, 26.0, 25.7), (16, 19.4, 23.2, 24.6), (48, 8.7, 22.8, 25.4)]

C_V1, C_V2, C_LOG = "#9e9e9e", "#3b6fb6", "#e08a2e"
W = 1000
X0 = 300               # bar origin
out: list[str] = []


def markup(s):
    """Escape ``s``; ``*...*`` becomes italic (species names)."""
    parts = escape(s).split("*")
    return "".join(f'<tspan font-style="italic">{p}</tspan>' if i % 2 else p
                   for i, p in enumerate(parts))


def text(x, y, s, size=12, anchor="start", color="#222", weight=None, italic=False):
    wt = f' font-weight="{weight}"' if weight else ""
    it = ' font-style="italic"' if italic else ""
    out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" '
               f'text-anchor="{anchor}"{wt}{it}>{markup(s)}</text>')


def rect(x, y, w, h, color):
    out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0.8):.1f}" height="{h}" fill="{color}"/>')


def line(x1, y1, x2, y2, color="#333", width=1.0, dash=False):
    da = ' stroke-dasharray="3 3"' if dash else ""
    out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
               f'stroke="{color}" stroke-width="{width}"{da}/>')


def fmt(minutes):
    return f"{minutes / 60:.1f} h" if minutes >= 90 else f"{minutes:.0f} min"


def axis(y_top, y_bot, max_min, step, scale, at):
    """Grid lines go in at index ``at`` (behind the bars), labels on top."""
    grid = []
    for t in range(0, max_min + 1, step):
        x = X0 + t * scale
        grid.append(f'<line x1="{x:.1f}" y1="{y_top:.1f}" x2="{x:.1f}" y2="{y_bot:.1f}" '
                    f'stroke="#e3e3e3" stroke-width="1"/>')
        text(x, y_bot + 15, f"{t}", size=11, anchor="middle", color="#555")
    out[at:at] = grid
    text(X0 + max_min * scale / 2, y_bot + 32, "wall time (min)", size=11.5,
         anchor="middle", color="#444")


def bar(y, h, logan, run, rng, color, scale, label_extra=""):
    rect(X0, y, logan * scale, h, C_LOG) if logan else None
    rect(X0 + logan * scale, y, run * scale, h, color)
    total = logan + run
    end = X0 + total * scale
    if rng:
        lo, hi = X0 + rng[0] * scale, X0 + rng[1] * scale
        line(lo, y + h / 2, hi, y + h / 2, "#222", 1.1)
        line(lo, y + 3, lo, y + h - 3, "#222", 1.1)
        line(hi, y + 3, hi, y + h - 3, "#222", 1.1)
        end = max(end, hi)
    text(end + 6, y + h / 2 + 4, fmt(total) + label_extra, size=11.5, weight="600")
    return end


# ---------------------------------------------------------------- panel A
y = 30
text(20, y, "A   Wall time per genome, 1000 batches, 48 threads", size=15, weight="700")
y += 14
scaleA = (W - X0 - 150) / 280
row_h, row_gap, group_gap = 18, 6, 18
y_axis_top = y + 8
at = len(out)
y += 16
for name, runs, rows in GENOMES:
    text(20, y + 12, name, size=12.5, weight="600", italic=True)
    text(20, y + 27, runs, size=11, color="#666")
    base = None
    for label, logan, run, rng, note in rows:
        text(X0 - 8, y + row_h / 2 + 4, label, size=11.5, anchor="end", color="#333")
        color = C_V1 if label.startswith("v1") else C_V2
        total = logan + run
        extra = ""
        if base and not label.startswith("v1"):
            extra = f"  ({base / total:.0f}× faster)"
        end = bar(y, row_h, logan, run, rng, color, scaleA, extra)
        if note:
            nx = end + 6 + 7.0 * len(fmt(total) + extra) + 10
            text(nx, y + row_h / 2 + 4, note, size=11, color="#666")
        if label.startswith("v1"):
            base = total
        y += row_h + row_gap
    y += group_gap
axis(y_axis_top, y - group_gap + 2, 280, 30, scaleA, at)
y += 30

# legend
lx = X0
for color, label in [(C_V1, "v1 behaviour: serial downloads, one final merge (new code, legacy flags)"),
                     (C_V2, "varus run (v2 defaults)"), (C_LOG, "varus logan pre-screen")]:
    rect(lx, y, 14, 12, color)
    text(lx + 20, y + 10.5, label, size=11.5, color="#333")
    y += 18
text(X0, y + 8, "Whiskers: range over 3 seeds. Mouse has no v1 run. S = VARUS score relative to A0 "
                "(mouse: intron sensitivity vs RefSeq).", size=11, color="#666")
y += 40

# ---------------------------------------------------------------- panel B
text(20, y, "B   Threads, *Chlorella sorokiniana* (seed 1)", size=15, weight="700")
y += 14
scaleB = (W - X0 - 150) / 110
y_axis_top = y + 8
at = len(out)
y += 16
for t, logan, run_l, run_n in THREADS:
    text(20, y + 22, f"--threads {t}", size=12.5, weight="600")
    text(X0 - 8, y + row_h / 2 + 4, "v2", size=11.5, anchor="end", color="#333")
    bar(y, row_h, 0, run_n, None, C_V2, scaleB)
    y += row_h + row_gap
    text(X0 - 8, y + row_h / 2 + 4, "v2 + Logan", size=11.5, anchor="end", color="#333")
    bar(y, row_h, logan, run_l, None, C_V2, scaleB)
    y += row_h + row_gap + 12
axis(y_axis_top, y - 10, 110, 10, scaleB, at)
y += 40
text(20, y, "Same results at every thread count. The read loop waits for downloads; "
            "the Logan pre-screen is limited by minimap2's CPU time.", size=11.5, color="#444")
H = int(y + 20)

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">\n'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>\n' + "\n".join(out) + "\n</svg>\n")
Path(__file__).with_name("runtime.svg").write_text(svg, encoding="utf-8")
