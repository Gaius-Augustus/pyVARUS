"""Draw docs/figures/threads.svg: total wall time with the Logan pre-screen by --threads.

    python docs/figures/threads.py      # writes threads.svg next to this file

Numbers are from docs/benchmark_logan.md, "Thread scaling" (Chlorella
sorokiniana GCA_025917655.1, 38.8 Mbp, 391 runs, 1000 batches, seed 1). Each bar is `varus logan` plus
`varus run` (including the final merge). Gray bars are the v1 loop on the
same genome: the old code in BRAKER4 production (2 threads, 4.27 h) and the
v1 behaviour of the new code (A0, 48 threads, 4.46 h).
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape

# (threads, (logan stage minutes, varus run minutes) or None, v1 minutes or None)
THREADS = [(2, None, 256.2), (4, (64.7, 37.1), None), (8, (35.6, 26.0), None),
           (16, (19.4, 23.2), None), (48, (8.7, 22.8), 267.6)]

C_RUN, C_LOG, C_V1 = "#3b6fb6", "#e08a2e", "#9e9e9e"
W, H = 660, 440
X0, Y0 = 80, 370       # plot origin (bottom left)
PW, PH = 540, 290      # plot width, height
Y_MAX, Y_STEP = 280, 40
out: list[str] = []


def markup(s):
    """Escape ``s``; ``*...*`` becomes italic (species names)."""
    parts = escape(s).split("*")
    return "".join(f'<tspan font-style="italic">{p}</tspan>' if i % 2 else p
                   for i, p in enumerate(parts))


def text(x, y, s, size=12, anchor="start", color="#222", weight=None, rotate=None):
    wt = f' font-weight="{weight}"' if weight else ""
    rot = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate else ""
    out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" '
               f'text-anchor="{anchor}"{wt}{rot}>{markup(s)}</text>')


def rect(x, y, w, h, color):
    out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{color}"/>')


def ypos(minutes):
    return Y0 - minutes / Y_MAX * PH


for t in range(0, Y_MAX + 1, Y_STEP):
    y = ypos(t)
    out.append(f'<line x1="{X0}" y1="{y:.1f}" x2="{X0 + PW}" y2="{y:.1f}" '
               f'stroke="#e3e3e3" stroke-width="1"/>')
    text(X0 - 8, y + 4, f"{t}", size=11.5, anchor="end", color="#555")
text(24, Y0 - PH / 2, "total wall time (min)", size=12.5, anchor="middle",
     color="#444", rotate=-90)

slot = PW / len(THREADS)
bw = slot * 0.36


def label(x, minutes):
    s = f"{minutes / 60:.1f} h" if minutes >= 120 else f"{minutes:.0f} min"
    text(x + bw / 2, ypos(minutes) - 7, s, size=12.5, anchor="middle", weight="600")


for i, (threads, v2, v1) in enumerate(THREADS):
    n = (v2 is not None) + (v1 is not None)
    x = X0 + slot * i + (slot - n * bw - (n - 1) * 4) / 2
    if v2:
        logan, run = v2
        rect(x, ypos(run), bw, Y0 - ypos(run), C_RUN)
        rect(x, ypos(run + logan), bw, ypos(run) - ypos(run + logan), C_LOG)
        label(x, logan + run)
        x += bw + 4
    if v1:
        rect(x, ypos(v1), bw, Y0 - ypos(v1), C_V1)
        label(x, v1)
    text(X0 + slot * (i + 0.5), Y0 + 18, f"{threads}", size=12.5, anchor="middle")
out.append(f'<line x1="{X0}" y1="{Y0}" x2="{X0 + PW}" y2="{Y0}" stroke="#333" stroke-width="1"/>')
text(X0 + PW / 2, Y0 + 40, "--threads", size=12.5, anchor="middle", color="#444")

text(20, 26, "1000 batches, *Chlorella sorokiniana* (38.8 Mbp), 391 runs",
     size=14, weight="700")
lx, ly = X0 + PW / 2 - 90, 50
for color, name in [(C_LOG, "v2: Logan pre-screen"), (C_RUN, "v2: online sampling"),
                    (C_V1, "v1: serial loop, no Logan")]:
    rect(lx, ly, 14, 12, color)
    text(lx + 20, ly + 10.5, name, size=12, color="#333")
    ly += 20

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">\n'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>\n' + "\n".join(out) + "\n</svg>\n")
Path(__file__).with_name("threads.svg").write_text(svg, encoding="utf-8")
