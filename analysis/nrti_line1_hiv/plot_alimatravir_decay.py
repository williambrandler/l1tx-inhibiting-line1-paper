"""One 11 mg alimatravir dose, intracellular triphosphate through day 38.

The peak is the phase 2 geometric-mean Cmax after 12 mg, scaled by 11/12
(Mayer et al., 2026, Table 3). The trough is the mean of the week 8, 12, 16,
and 24 arithmetic means in the same study (Table 4), scaled the same way.
Those troughs are taken about 28 days after a dose. The half-life is the
single exponential that joins the peak, at 1 day, to that trough. Day 38 is
that exponential continued past the last measured trough.

The HIV line is the published 0.03 pmol/10^6 threshold for a ≥1 log10 drop
in viral RNA. Alimatravir has no published LINE-1 IC50. The Figure 2A fit
predicts 1.3 nM from the 0.21 nM HIV IC50, and the intracellular amount at
that IC50 is 9.2 fmol/10^6 times 1.3/0.21. That concentration stops half a
full reporter cassette. The LINE-1 lines scale it to the triphosphate at
which 90% of 100-nt products, and 90% of 1,000-nt products, have stopped.
The first 40 nt are the polyA tail and contribute no adenine-analogue site.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT.parents[1] / "paper" / "line1-inhibition"
OUT = OUT_DIR / "alimatravir_11mg_decay.svg"
PNG = OUT_DIR / "alimatravir_11mg_decay.png"

PEAK = 1.48 * 11 / 12
TROUGH_12 = (0.0815 + 0.0715 + 0.0882 + 0.0735) / 4
TROUGH = TROUGH_12 * 11 / 12
HALF_LIFE_H = 27 * 24 * math.log(2) / -math.log(TROUGH / PEAK)
HIV_IC50_NM = 0.21
LINE1_IC50_NM = 1.3
IC50 = 0.0092 * LINE1_IC50_NM / HIV_IC50_NM
HIV90 = 0.03
N_T = 497


def line1_90(length_nt: float) -> float:
    """90% stop for a product whose first 40 nt are polyA."""
    n_sites = N_T * (length_nt - 40) / 1803
    return (-math.log(0.1) * N_T / (n_sites * math.log(2))) * IC50


LINE1_90 = line1_90(100)
LINE1_1000 = line1_90(1000)

INK = "#000000"


def concentration(day: float) -> float:
    if day <= 1.0:
        return PEAK
    return PEAK * 0.5 ** ((day - 1.0) * 24.0 / HALF_LIFE_H)


def day_at(level: float) -> float:
    return 1.0 + HALF_LIFE_H / 24.0 * math.log(level / PEAK) / math.log(0.5)


def _map(value: float, src: tuple[float, float], dst: tuple[float, float]) -> float:
    low, high = src
    y0, y1 = dst
    return y0 + (value - low) / (high - low) * (y1 - y0)


WIDTH = 960
HEIGHT = 640


def build_svg() -> str:
    width, height = WIDTH, HEIGHT
    left, right, top, bottom = 168, 44, 40, 100
    plot_right = width - right
    plot_bottom = height - bottom
    x_max = 38.0
    y_min, y_max = math.log10(0.008), math.log10(20.0)

    def px(day: float) -> float:
        return _map(day, (0.0, x_max), (left, plot_right))

    def py(level: float) -> float:
        return _map(math.log10(level), (y_min, y_max), (plot_bottom, top))

    steps = 380
    series = []
    for index in range(steps + 1):
        day = 1.0 + 37.0 * index / steps
        series.append((px(day), py(concentration(day))))
    points = " ".join(f"{x:.1f},{y:.1f}" for x, y in series)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif" fill="{INK}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{plot_bottom}" stroke="{INK}" stroke-width="1.2"/>',
        f'<line x1="{left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="{INK}" stroke-width="1.2"/>',
    ]

    for tick in (0.01, 0.1, 1, 10):
        y = py(tick)
        parts.append(
            f'<line x1="{left - 8}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="{INK}"/>'
        )
        parts.append(
            f'<text x="{left - 14}" y="{y + 8:.1f}" text-anchor="end" font-size="24">{tick:g}</text>'
        )

    guides = (
        (LINE1_90, "6 4", "LINE-1, 90% inhibition of 100-nt products", 6, -12, "start"),
        (LINE1_1000, "6 4", "LINE-1, 1000-nt", 38, -16, "end"),
        (HIV90, "1.2 3.2", "HIV, 90% inhibition of 9.7-kb products", 2, -12, "start"),
    )
    for level, dash, label, label_day, label_dy, anchor in guides:
        y = py(level)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" '
            f'stroke="{INK}" stroke-width="1.4" stroke-dasharray="{dash}"/>'
        )
        parts.append(
            f'<text x="{px(label_day):.1f}" y="{y + label_dy:.1f}" text-anchor="{anchor}" '
            f'font-size="24" font-weight="bold">{label}</text>'
        )

    parts.append(
        f'<polyline fill="none" stroke="{INK}" stroke-width="2.2" '
        f'stroke-linejoin="round" stroke-linecap="round" points="{points}"/>'
    )

    markers = (
        (28.0, TROUGH, "measured trough", 14, -18, "start"),
        (38.0, concentration(38.0), "", 0, 0, "start"),
    )
    for day, level, label, dx, dy, anchor in markers:
        x, y = px(day), py(level)
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{INK}" stroke="#ffffff" stroke-width="1.2"/>'
        )
        if label:
            parts.append(
                f'<text x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}" font-size="22">{label}</text>'
            )

    for tick in (0, 7, 14, 21, 28, 38):
        x = px(tick)
        parts.append(
            f'<line x1="{x:.1f}" y1="{plot_bottom}" x2="{x:.1f}" y2="{plot_bottom + 8}" stroke="{INK}"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{plot_bottom + 36}" text-anchor="middle" font-size="24">{tick}</text>'
        )

    parts.append(
        f'<text x="{(left + plot_right) / 2:.1f}" y="{plot_bottom + 74}" text-anchor="middle" '
        f'font-size="26" font-weight="bold">Days after dose</text>'
    )
    ylabel_y = (top + plot_bottom) / 2
    parts.append(
        f'<text x="28" y="{ylabel_y:.1f}" font-size="26" font-weight="bold" text-anchor="middle" '
        f'transform="rotate(-90 28 {ylabel_y:.1f})">Alimatravir triphosphate (pmol/10\u2076 cells)</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def write_png(svg_path: Path, png_path: Path, width: int, height: int) -> None:
    """Rasterize the SVG at 3× so a ~7 inch figure is near 300 dpi."""
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
    svg = build_svg()
    OUT.write_text(svg)
    write_png(OUT, PNG, WIDTH, HEIGHT)
    print(f"peak {PEAK:.3f}")
    print(f"trough {TROUGH:.3f}")
    print(f"half-life {HALF_LIFE_H:.1f} h")
    for name, level in (("100-nt", LINE1_90), ("1000-nt", LINE1_1000)):
        if level >= PEAK:
            print(f"line1 {name} {level:.3f} above the peak")
        else:
            print(f"line1 {name} {level:.3f} at day {day_at(level):.1f}")
    print(f"hiv 90 {HIV90:.3f} at day {day_at(HIV90):.1f}")
    print(f"day 38 {concentration(38):.3f}")
    print(OUT)
    print(PNG)


if __name__ == "__main__":
    main()
