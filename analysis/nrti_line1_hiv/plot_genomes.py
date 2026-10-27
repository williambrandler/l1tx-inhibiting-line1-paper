"""Two-panel figure. A, the LINE-1 and HIV-1 genomes. B, what an insertion looks like.

Panel A is drawn from HXB2 coordinates for HIV-1 (K03455) and L1.3 coordinates
for LINE-1 (L19088), both at one scale.

Panel B is the structure left in genomic DNA after target-primed reverse
transcription. The classical product keeps the 3' end of the element, a polyA
tail, and a target-site duplication from the staggered endonuclease cut; the 5'
end is lost, leaving a median insert near 1 kb of the 6 kb element (Ostertag and
Kazazian, 2001; Lander et al., 2001). PolyA tails on young insertions run about
20-100 bp (Grandi et al., 2013). The non-classical product is truncated at
both ends, carries no polyA tail, and has no target-site duplication
(Morrish et al., 2002; Sen et al., 2007).

The non-classical fragment is drawn at the size these two datasets report, but
unlabelled, because neither is a clean estimate. The 23 PacBio-validated
3'-truncated brain insertions in Jacob-Hirsch et al. (2018), Table S4, run
235-3 856 bp with a median near 650 bp, and their 3' ends sit mostly between
L1Hs positions 4 500 and 5 500; the 21 endonuclease-independent loci in
Sen et al. (2007) run 34-4 410 bp.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT.parents[1] / "paper" / "line1-inhibition"
OUT = OUT_DIR / "genomes.svg"
PNG = OUT_DIR / "genomes.png"

WIDTH = 980
HEIGHT = 616

INK = "#1a1a1a"
GREY = "#8a8a8a"
DNA = "#6e6e6e"
UTR = "#d9d9d9"
ORF1 = "#f4c7a8"
ORF2 = "#f0a0a0"
POLYA = "#9fc7e8"
TSD = "#3f3f3f"

FONT = "Helvetica, Arial, sans-serif"

PX_PER_KB_A = 88.0
PX_PER_KB_B = 118.0

# LINE-1, L1.3 numbering, kb.
L1_FEATURES = [
    (0.000, 0.907, "5' UTR", UTR, "in"),
    (0.907, 1.918, "ORF1", ORF1, "in"),
    (1.918, 5.796, "ORF2, EN / RT", ORF2, "in"),
    (5.796, 6.000, "3'", UTR, "above"),
]

# HIV-1, HXB2 numbering, kb. Order is draw order: wide genes first.
HIV_FEATURES = [
    (0.001, 0.634, "LTR", UTR, "in"),
    (0.790, 2.292, "gag", ORF1, "in"),
    (2.085, 5.096, "pol (RT)", ORF2, "in"),
    (5.041, 5.619, "vif", "#e7d48b", "in"),
    (5.559, 5.850, "vpr", "#c9c4e0", "above1"),
    (5.831, 6.045, "tat", "#9eb0d8", "above2"),
    (5.970, 6.045, "rev", "#f3c97a", "above1"),
    (6.062, 6.310, "vpu", "#f6e27a", "above2"),
    (6.225, 8.795, "env", "#b7d3c6", "in"),
    (8.379, 8.469, "tat", "#9eb0d8", "above1"),
    (8.379, 8.653, "rev", "#f3c97a", "above2"),
    (8.797, 9.417, "nef", "#d7c2de", "in"),
    (9.086, 9.719, "LTR", UTR, "in"),
]

# nef is half covered by the 3' LTR, so its label sits over the visible part.
LABEL_KB = {"nef": 8.94}

HIV_LEN = 9.719
L1_LEN = 6.000

# Panel B inserts, kb of the element that survives.
CLASSICAL = (5.000, 6.000)
NONCLASSICAL = (4.350, 5.000)


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def jag_path(x0: float, x1: float, y0: float, y1: float, left: bool, right: bool) -> str:
    """Box outline, with a torn zigzag on the truncated ends."""
    teeth, amp = 5, 4.0
    out = [f"M {x0:.1f} {y0:.1f}", f"L {x1:.1f} {y0:.1f}"]
    if right:
        for i in range(1, teeth + 1):
            y = y0 + (y1 - y0) * i / teeth
            dx = amp if i % 2 else -amp
            out.append(f"L {x1 + (0 if i == teeth else dx):.1f} {y:.1f}")
    else:
        out.append(f"L {x1:.1f} {y1:.1f}")
    out.append(f"L {x0:.1f} {y1:.1f}")
    if left:
        for i in range(1, teeth + 1):
            y = y1 - (y1 - y0) * i / teeth
            dx = amp if i % 2 else -amp
            out.append(f"L {x0 + (0 if i == teeth else dx):.1f} {y:.1f}")
    out.append("Z")
    return " ".join(out)


def gene_track(
    features: list[tuple[float, float, str, str, str]],
    x0: float,
    y: float,
    height: float,
    px_per_kb: float,
    total_kb: float,
) -> list[str]:
    parts = [
        f'<rect x="{x0:.1f}" y="{y:.1f}" width="{total_kb * px_per_kb:.1f}" '
        f'height="{height:.1f}" fill="#fbfbfb" stroke="{INK}" stroke-width="1"/>'
    ]
    for start, end, _label, fill, _place in features:
        bx = x0 + start * px_per_kb
        bw = (end - start) * px_per_kb
        parts.append(
            f'<rect x="{bx:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{height:.1f}" '
            f'fill="{fill}" fill-opacity="0.92" stroke="{INK}" stroke-width="1"/>'
        )
    for start, end, label, _fill, place in features:
        mid = x0 + LABEL_KB.get(label, (start + end) / 2) * px_per_kb
        if place == "in":
            parts.append(
                f'<text x="{mid:.1f}" y="{y + height / 2 + 4:.1f}" text-anchor="middle" '
                f'font-size="12">{esc(label)}</text>'
            )
        else:
            lift = 10 if place in ("above", "above1") else 25
            parts.append(
                f'<line x1="{mid:.1f}" y1="{y - 2:.1f}" x2="{mid:.1f}" y2="{y - lift + 1:.1f}" '
                f'stroke="{GREY}" stroke-width="0.8"/>'
            )
            parts.append(
                f'<text x="{mid:.1f}" y="{y - lift - 2:.1f}" text-anchor="middle" '
                f'font-size="11">{esc(label)}</text>'
            )
    return parts


def retained(x0: float, y: float, height: float, span: tuple[float, float]) -> list[str]:
    """The piece of the element that is actually in the genome, in element colours."""
    start, end = span
    parts = []
    for f_start, f_end, _label, fill, _place in L1_FEATURES:
        lo, hi = max(f_start, start), min(f_end, end)
        if hi <= lo:
            continue
        bx = x0 + lo * PX_PER_KB_B
        bw = (hi - lo) * PX_PER_KB_B
        parts.append(
            f'<rect x="{bx:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{height:.1f}" '
            f'fill="{fill}" fill-opacity="0.92" stroke="none"/>'
        )
    bx0 = x0 + start * PX_PER_KB_B
    bx1 = x0 + end * PX_PER_KB_B
    left_jag = start > 0.001
    right_jag = end < L1_LEN - 0.001
    parts.append(
        f'<path d="{jag_path(bx0, bx1, y, y + height, left_jag, right_jag)}" '
        f'fill="none" stroke="{INK}" stroke-width="1.2"/>'
    )
    return parts


def dna(x0: float, x1: float, y: float) -> str:
    return (
        f'<rect x="{x0:.1f}" y="{y - 3:.1f}" width="{x1 - x0:.1f}" height="6" '
        f'fill="{DNA}" stroke="none"/>'
    )


def tsd_box(x: float, y: float, height: float) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="11" height="{height:.1f}" '
        f'fill="{TSD}" stroke="{INK}" stroke-width="1"/>'
    )


def build_svg() -> str:
    p: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="{FONT}" fill="{INK}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]

    # Panel A.
    ax = 48.0
    p.append(f'<text x="18" y="30" font-size="19" font-weight="700">A</text>')
    p.append(
        f'<text x="48" y="30" font-size="15" font-weight="700">'
        f'Structure of LINE-1 and HIV genomes</text>'
    )
    p.append(f'<text x="48" y="58" font-size="13" font-weight="700">LINE-1, 6.0 kb, two genes</text>')
    p += gene_track(L1_FEATURES, ax, 68.0, 26.0, PX_PER_KB_A, L1_LEN)
    p.append(
        f'<text x="48" y="136" font-size="13" font-weight="700">'
        f'HIV-1, 9.7 kb, nine genes; tat and rev are each two exons</text>'
    )
    p += gene_track(HIV_FEATURES, ax, 182.0, 26.0, PX_PER_KB_A, HIV_LEN)

    # Panel B.
    p.append(f'<text x="18" y="272" font-size="19" font-weight="700">B</text>')
    p.append(
        f'<text x="48" y="272" font-size="15" font-weight="700">'
        f'Structure of LINE-1 insertions</text>'
    )

    bx = 110.0
    ref_y = 300.0
    h = 24.0
    p.append(
        f'<text x="48" y="{ref_y - 10:.1f}" font-size="12" fill="{GREY}">'
        f'full-length element, 6.0 kb</text>'
    )
    for f_start, f_end, _label, fill, _place in L1_FEATURES:
        p.append(
            f'<rect x="{bx + f_start * PX_PER_KB_B:.1f}" y="{ref_y:.1f}" '
            f'width="{(f_end - f_start) * PX_PER_KB_B:.1f}" height="{h:.1f}" '
            f'fill="{fill}" fill-opacity="0.30" stroke="{GREY}" stroke-width="0.8" '
            f'stroke-dasharray="3 2"/>'
        )
    for kb in range(7):
        tick = bx + kb * PX_PER_KB_B
        p.append(
            f'<line x1="{tick:.1f}" y1="{ref_y + h:.1f}" x2="{tick:.1f}" '
            f'y2="{ref_y + h + 5:.1f}" stroke="{GREY}"/>'
        )
        p.append(
            f'<text x="{tick:.1f}" y="{ref_y + h + 17:.1f}" text-anchor="middle" '
            f'font-size="10" fill="{GREY}">{kb}</text>'
        )
    p.append(
        f'<text x="{bx + 6 * PX_PER_KB_B + 16:.1f}" y="{ref_y + h + 17:.1f}" '
        f'font-size="10" fill="{GREY}">kb</text>'
    )

    rows = [
        ("Classical, endonuclease-dependent", CLASSICAL, True, 428.0),
        ("Non-classical, endonuclease-independent", NONCLASSICAL, False, 552.0),
    ]
    for title, span, classical, row_y in rows:
        start, end = span
        x_ins0 = bx + start * PX_PER_KB_B
        x_ins1 = bx + end * PX_PER_KB_B
        mid_y = row_y + h / 2
        p.append(f'<text x="48" y="{row_y - 34:.1f}" font-size="13" font-weight="700">{esc(title)}</text>')
        p.append(
            f'<text x="52" y="{mid_y - 10:.1f}" font-size="11" fill="{DNA}">genomic DNA</text>'
        )
        # Dotted drops from the reference element to the ends that survive.
        for x_drop in (x_ins0, x_ins1):
            p.append(
                f'<line x1="{x_drop:.1f}" y1="{ref_y + h + 22:.1f}" x2="{x_drop:.1f}" '
                f'y2="{row_y - 6:.1f}" stroke="{GREY}" stroke-width="0.8" stroke-dasharray="2 3"/>'
            )

        if classical:
            tail_w = 44.0
            x_tail1 = x_ins1 + tail_w
            x_tsd0 = x_ins0 - 17
            p.append(dna(48.0, x_tsd0, mid_y))
            p.append(dna(x_tail1 + 17, WIDTH - 40, mid_y))
            p.append(tsd_box(x_tsd0, row_y, h))
            p.append(tsd_box(x_tail1 + 6, row_y, h))
            p += retained(bx, row_y, h, span)
            p.append(
                f'<rect x="{x_ins1:.1f}" y="{row_y:.1f}" width="{tail_w:.1f}" height="{h:.1f}" '
                f'fill="{POLYA}" fill-opacity="0.92" stroke="{INK}" stroke-width="1.2"/>'
            )
            p.append(
                f'<text x="{x_ins1 + tail_w / 2:.1f}" y="{mid_y + 4:.1f}" text-anchor="middle" '
                f'font-size="11">(A)n</text>'
            )
            for x_t in (x_tsd0 + 5.5, x_tail1 + 11.5):
                p.append(
                    f'<text x="{x_t:.1f}" y="{row_y + h + 16:.1f}" text-anchor="middle" '
                    f'font-size="11">TSD</text>'
                )
            p.append(
                f'<text x="{(x_ins0 + x_ins1) / 2:.1f}" y="{row_y - 12:.1f}" text-anchor="middle" '
                f'font-size="11">median insert \u2248 1 kb</text>'
            )
            p.append(
                f'<text x="{x_ins0 - 6:.1f}" y="{row_y - 30:.1f}" text-anchor="end" '
                f'font-size="11" fill="{GREY}">5\u2032 truncated</text>'
            )
            p.append(
                f'<text x="{WIDTH - 44}" y="{row_y + h + 32:.1f}" '
                f'text-anchor="end" font-size="11" fill="{GREY}">'
                f'polyA tail (20\u2013100 bp), target site duplication (7\u201320 bp) at both '
                f'junctions</text>'
            )
        else:
            p.append(dna(48.0, x_ins0, mid_y))
            p.append(dna(x_ins1, WIDTH - 40, mid_y))
            p += retained(bx, row_y, h, span)
            p.append(
                f'<text x="{(x_ins0 + x_ins1) / 2:.1f}" y="{row_y - 12:.1f}" text-anchor="middle" '
                f'font-size="11">internal fragment, both ends lost</text>'
            )
            p.append(
                f'<text x="{x_ins0 - 16:.1f}" y="{row_y - 30:.1f}" text-anchor="end" '
                f'font-size="11" fill="{GREY}">5\u2032 truncated</text>'
            )
            p.append(
                f'<text x="{x_ins1 + 16:.1f}" y="{row_y - 30:.1f}" font-size="11" '
                f'fill="{GREY}">3\u2032 truncated</text>'
            )
            p.append(
                f'<text x="{(x_ins0 + x_ins1) / 2:.1f}" y="{row_y + h + 16:.1f}" '
                f'text-anchor="middle" font-size="11">no polyA tail, no target site duplication</text>'
            )

    p.append("</svg>")
    return "\n".join(p)


def write_png(svg_path: Path, png_path: Path, width: int, height: int) -> None:
    """Rasterize the SVG at 3x so a ~7 inch figure is near 300 dpi."""
    scale = 3
    pixel_w, pixel_h = width * scale, height * scale
    svg = svg_path.read_text()
    svg = svg.replace(f'width="{width}"', f'width="{pixel_w}"', 1)
    svg = svg.replace(f'height="{height}"', f'height="{pixel_h}"', 1)
    html_path = png_path.with_suffix(".html")
    html_path.write_text(
        "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
        f"<body style='margin:0;background:#ffffff;width:{pixel_w}px;height:{pixel_h}px'>"
        f"{svg}</body></html>"
    )
    chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    subprocess.run(
        [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=1",
            f"--window-size={pixel_w},{pixel_h}",
            f"--screenshot={png_path}",
            html_path.as_uri(),
        ],
        check=True,
        capture_output=True,
    )
    html_path.unlink(missing_ok=True)


def main() -> None:
    OUT.write_text(build_svg())
    write_png(OUT, PNG, WIDTH, HEIGHT)
    print(OUT)
    print(PNG)


if __name__ == "__main__":
    main()
