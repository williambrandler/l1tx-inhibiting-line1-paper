"""LINE-1 versus HIV potency. One point per nucleoside.

Cell-assay IC50s stay on the nucleoside panels. Enzyme Ki and IC50 values
stay on the triphosphate panels. Where several assays reported the same
nucleoside, the point is their geometric mean.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "potency_data.tsv"
OUT = ROOT / "line1_hiv_correlation.html"
PER_KB_OUT = ROOT / "line1_hiv_per_kb.html"
PNG_OUT = ROOT.parents[1] / "paper" / "line1-inhibition" / "line1_hiv_correlation.png"

# HXB2 (K03455) is the HIV genome the cell assays have to finish.
# The LINE-1 number is the pYX017 HeLa assay, which scores the cassette in
# the 3' UTR: 150 bp of the L1.3 3' end, a 40-base polyA, and the 1,653 bp
# firefly luciferase coding sequence. Counts are A, T, G, C.
HIV_BASES = {"A": 3411, "T": 2163, "G": 2373, "C": 1772}
LINE1_BASES = {"A": 591, "T": 497, "G": 394, "C": 361}
# An analog is inserted where that natural base would go, so the template
# base is its complement. Adenine opposite T, thymine opposite A.
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
NUCLEOSIDE_BASE = {
    "EFdA": "A",
    "2-Cl-EdA": "A",
    "EdA": "A",
    "EtAdA": "A",
    "EdDAP": "A",
    "EdG": "G",
    "EdC": "C",
    "4'-E-FdC": "C",
    "TAF": "A",
    "TDF": "A",
    "TFV": "A",
    "ABC": "G",
    "AZT": "T",
    "d4T": "T",
    "STP": "T",
    "4'-Ed4T": "T",
    "4'-E-T": "T",
    "3TC": "C",
    "FTC": "C",
    "L-d4FC": "C",
    "ATC": "C",
}

# 0.4 pL per PBMC. 1 pmol/10^6 cells = 2.5 µM = 2500 nM.
PMOL_TO_NM = 2500.0
# Relative potency is 1 at 100 µM on both axes.
REFERENCE_NM = 100_000.0

NUCLEOSIDE_TYPES = {"nucleoside", "nucleoside_jones", "nucleoside_baldwin"}
TRIPHOSPHATE_TYPES = {
    "triphosphate",
    "triphosphate_baldwin_htrf",
    "triphosphate_baldwin_gel",
    "triphosphate_baldwin_fl",
}


def load() -> list[dict]:
    with DATA.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key in (
            "line1_nM",
            "hiv_nM",
            "phenosense_hiv_nM",
            "dose_mg",
            "cmax_nM",
            "pbmc_pmol_per_1e6",
        ):
            text = (row.get(key) or "").strip()
            row[key] = float(text) if text else None
    return rows


def geomean(values: list[float]) -> float:
    cleaned = [value for value in values if value and value > 0]
    return float(np.exp(np.mean(np.log(cleaned))))


def log_fit(x: np.ndarray, y: np.ndarray) -> dict:
    lx, ly = np.log10(x), np.log10(y)
    slope, intercept, r_pearson, p_pearson, _ = stats.linregress(lx, ly)
    rho, p_rho = stats.spearmanr(lx, ly)
    return {
        "slope": slope,
        "intercept": intercept,
        "r": r_pearson,
        "p": p_pearson,
        "rho": rho,
        "p_rho": p_rho,
        "n": len(x),
    }


def decade_ticks(values: np.ndarray) -> tuple[list[float], list[str], list[float]]:
    """Ticks at 1, 10, 100, … labeled 0, 1, 10, 100, …"""
    hi = int(np.floor(np.log10(values.max() * 1.35)))
    vals = [float(10**power) for power in range(0, hi + 1)]
    text = ["0" if value == 1 else str(int(value / 10)) for value in vals]
    span = [np.log10(0.65), np.log10(values.max() * 1.7)]
    return vals, text, span


def fit_line(fit: dict, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grid = np.logspace(np.log10(x.min() / 1.5), np.log10(x.max() * 1.5), 40)
    y = 10 ** (fit["intercept"] + fit["slope"] * np.log10(grid))
    return grid, y


def collapse(rows: list[dict], types: set[str]) -> list[dict]:
    """One point per short name. LINE-1 is the geometric mean across assays."""
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        if row["record_type"] not in types or not row["line1_nM"] or not row["hiv_nM"]:
            continue
        grouped.setdefault(row["short"], []).append(row)
    points = []
    for short, members in grouped.items():
        primary = next(
            (row for row in members if row["record_type"] in {"nucleoside", "triphosphate"}),
            members[0],
        )
        point = dict(primary)
        point["short"] = short
        point["line1_nM"] = geomean([row["line1_nM"] for row in members])
        point["n_assays"] = len(members)
        point["line1_values"] = [row["line1_nM"] for row in members]
        if "greater than" in (primary.get("line1_assay") or "").lower():
            point["line1_censor"] = "gt"
        points.append(point)
    return points


def hover(point: dict, extra: str = "") -> str:
    values = ", ".join(f"{value:g}" for value in point["line1_values"])
    if point.get("line1_censor") == "gt":
        line1_line = f"LINE-1 &gt;{point['line1_nM']:.3g} nM, plotted at the bin edge"
    else:
        line1_line = f"LINE-1 {point['line1_nM']:.3g} nM, geometric mean of {point['n_assays']}"
    parts = [
        f"<b>{point['name']}</b>",
        line1_line,
        values + " nM",
        f"HIV {point['hiv_nM']:g} nM",
    ]
    if extra:
        parts.append(extra)
    return "<br>".join(parts)


def nucleoside_points() -> list[dict]:
    nucleosides = collapse(load(), NUCLEOSIDE_TYPES)
    nucleosides.sort(key=lambda point: point["hiv_nM"])
    return nucleosides


def abbreviations(nucleosides: list[dict]) -> list[tuple[str, str]]:
    ordered = sorted(nucleosides, key=lambda point: point["short"].lower())
    return [(point["short"], point["name"]) for point in ordered]


def per_kb_ic50(ic50_nM: float, bases: dict[str, int], nucleoside_base: str) -> float:
    """IC50 restated per 1,000 cognate template bases.

    A longer template, or one richer in the complementary base, makes the
    same per-site effect look like a lower IC50. Dividing that advantage
    out leaves the potency of 1,000 stops.
    """
    template = COMPLEMENT[nucleoside_base]
    return ic50_nM * bases[template] / 1000.0


def per_kb_points(nucleosides: list[dict]) -> list[dict]:
    adjusted = []
    for point in nucleosides:
        base = NUCLEOSIDE_BASE.get(point["short"])
        if base is None:
            raise KeyError(f"no base assigned for {point['short']}")
        copy = dict(point)
        copy["nucleoside_base"] = base
        copy["template_base"] = COMPLEMENT[base]
        copy["hiv_nM_raw"] = point["hiv_nM"]
        copy["line1_nM_raw"] = point["line1_nM"]
        copy["hiv_nM"] = per_kb_ic50(point["hiv_nM"], HIV_BASES, base)
        copy["line1_nM"] = per_kb_ic50(point["line1_nM"], LINE1_BASES, base)
        adjusted.append(copy)
    return adjusted


def per_kb_hover(point: dict) -> str:
    template = point["template_base"]
    hiv_sites = HIV_BASES[template]
    line1_sites = LINE1_BASES[template]
    line1_raw = f"LINE-1 {point['line1_nM_raw']:.3g} nM"
    if point.get("line1_censor") == "gt":
        line1_raw = f"LINE-1 &gt;{point['line1_nM_raw']:.3g} nM, bin edge"
    return "<br>".join(
        [
            f"<b>{point['name']}</b>",
            f"Base {point['nucleoside_base']}, incorporated opposite {template}",
            f"HIV {point['hiv_nM_raw']:.3g} nM over {sum(HIV_BASES.values())} bp, {hiv_sites} {template} sites",
            f"{line1_raw} over {sum(LINE1_BASES.values())} bp, {line1_sites} {template} sites",
            f"Per 1,000 {template} sites: HIV {point['hiv_nM']:.3g} nM, LINE-1 {point['line1_nM']:.3g} nM",
        ]
    )


def scatter_figure(
    nucleosides: list[dict],
    *,
    width: int,
    height: int,
    text_size: int = 12,
    marker_size: int = 11,
    axis_size: int = 13,
    tick_size: int = 12,
    color: str = "#feb24c",
    margin: dict | None = None,
    x_title: str = "HIV relative potency",
    y_title: str = "LINE-1 relative potency",
    hover_text=None,
) -> go.Figure:
    hiv = REFERENCE_NM / np.array([point["hiv_nM"] for point in nucleosides])
    line1 = REFERENCE_NM / np.array([point["line1_nM"] for point in nucleosides])
    fitted_idx = [i for i, point in enumerate(nucleosides) if point.get("line1_censor") != "gt"]
    censored_idx = [i for i, point in enumerate(nucleosides) if point.get("line1_censor") == "gt"]
    labels = hover_text or (lambda point: hover(point))
    fit_n = log_fit(hiv[fitted_idx], line1[fitted_idx])
    grid, fitted = fit_line(fit_n, hiv[fitted_idx])
    x_ticks, x_ticktext, x_range = decade_ticks(hiv)
    y_ticks, y_ticktext, y_range = decade_ticks(line1)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=grid,
            y=fitted,
            mode="lines",
            line=dict(color="#000000", width=1.5),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    def add_points(indexes: list[int], marker: dict) -> None:
        if not indexes:
            return
        chosen = [nucleosides[i] for i in indexes]
        fig.add_trace(
            go.Scatter(
                x=hiv[indexes],
                y=line1[indexes],
                mode="markers+text",
                text=[point["short"] for point in chosen],
                textposition="top right",
                textfont=dict(size=text_size, color="#000000", family="Inter, sans-serif"),
                cliponaxis=False,
                marker=marker,
                hovertext=[labels(point) for point in chosen],
                hoverinfo="text",
                showlegend=False,
            )
        )

    add_points(fitted_idx, dict(size=marker_size, color=color))
    add_points(censored_idx, dict(size=marker_size, color=color))
    fig.update_xaxes(
        type="log",
        title_text=x_title,
        range=x_range,
        tickmode="array",
        tickvals=x_ticks,
        ticktext=x_ticktext,
        minor=dict(ticks="", showgrid=False),
        title_font=dict(family="Inter, sans-serif", size=axis_size, color="#000000"),
        tickfont=dict(family="Inter, sans-serif", size=tick_size, color="#000000"),
        linecolor="#000000",
        tickcolor="#000000",
        automargin=False,
    )
    fig.update_yaxes(
        type="log",
        title_text=y_title,
        range=y_range,
        tickmode="array",
        tickvals=y_ticks,
        ticktext=y_ticktext,
        minor=dict(ticks="", showgrid=False),
        title_font=dict(family="Inter, sans-serif", size=axis_size, color="#000000"),
        tickfont=dict(family="Inter, sans-serif", size=tick_size, color="#000000"),
        linecolor="#000000",
        tickcolor="#000000",
        automargin=False,
    )
    fig.update_layout(
        height=height,
        width=width,
        template="plotly_white",
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font=dict(family="Helvetica, Arial, sans-serif", size=tick_size, color="#000000"),
        margin=margin or dict(t=24, b=50, l=70, r=40),
    )
    fig.add_annotation(
        x=0.03,
        y=0.97,
        xref="x domain",
        yref="y domain",
        text=f"r = {fit_n['r']:.2f}",
        showarrow=False,
        xanchor="left",
        yanchor="top",
        font=dict(family="Helvetica, Arial, sans-serif", size=axis_size, color="#000000"),
    )
    return fig


def write_per_kb() -> None:
    points = per_kb_points(nucleoside_points())
    fig = scatter_figure(
        points,
        width=760,
        height=680,
        x_title="HIV relative potency per 1 kb",
        y_title="LINE-1 relative potency per 1 kb",
        hover_text=per_kb_hover,
    )
    names = "".join(
        f"<div><span>{point['short']}</span>: {point['name']} ({point['nucleoside_base']})</div>"
        for point in sorted(points, key=lambda item: item["short"].lower())
    )
    body = fig.to_html(full_html=False, include_plotlyjs="cdn", default_width="760px")
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HIV and LINE-1 potency per 1 kb</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Space+Grotesk:wght@500;600&display=swap" rel="stylesheet">
<style>
  body {{ font-family: Inter, sans-serif; margin: 24px auto; max-width: 1180px; color: #000000; background: #ffffff; }}
  h1 {{ font-family: "Space Grotesk", sans-serif; font-weight: 500; color: #000000; }}
  p {{ max-width: 760px; line-height: 1.5; }}
  .row {{ display: flex; align-items: flex-start; gap: 36px; }}
  aside {{ color: #000000; font-size: 15px; line-height: 1.7; padding-top: 72px; }}
  aside span {{ font-weight: 600; }}
</style>
</head>
<body>
<h1>HIV &amp; LINE-1 potency per 1 kb</h1>
<p>Same nucleosides and the same relative-potency scale as the dashboard scatter. Each IC50 is restated per 1,000 template bases. HIV is HXB2, 9,719 bp. LINE-1 is the scored pYX017 stretch, 1,843 bp.</p>
<div class="row">
{body}
<aside>
{names}
</aside>
</div>
<script>
document.fonts.ready.then(function () {{
  document.querySelectorAll(".js-plotly-plot").forEach(function (el) {{
    Plotly.relayout(el, {{}});
  }});
}});
</script>
</body>
</html>
"""
    PER_KB_OUT.write_text(html)
    fitted = [point for point in points if point.get("line1_censor") != "gt"]
    hiv = REFERENCE_NM / np.array([point["hiv_nM"] for point in fitted])
    line1 = REFERENCE_NM / np.array([point["line1_nM"] for point in fitted])
    fit_n = log_fit(hiv, line1)
    print(f"wrote {PER_KB_OUT}")
    print(
        f"per kb n={fit_n['n']} slope={fit_n['slope']:.3f} "
        f"rho={fit_n['rho']:.3f} p={fit_n['p_rho']:.4g}"
    )
    for point in sorted(points, key=lambda item: item["line1_nM"] / item["hiv_nM"]):
        ratio = point["line1_nM"] / point["hiv_nM"]
        print(
            f"  {point['short']:8} {point['nucleoside_base']} "
            f"LINE-1/HIV per kb {ratio:.2f}"
        )


def write_png(nucleosides: list[dict]) -> None:
    """Static figure for the paper. Plotly JS is inlined so Chrome needs no network."""
    # Scatter is 80% of the previous paper figure. The abbreviation key
    # beside it keeps the same type size.
    fig = scatter_figure(
        nucleosides,
        width=608,
        height=544,
        text_size=13,
        marker_size=10,
        axis_size=17,
        tick_size=15,
        margin=dict(t=22, b=58, l=86, r=51),
    )
    legend = "".join(
        f"<div><span>{short}</span> = {name}</div>"
        for short, name in abbreviations(nucleosides)
    )
    body = fig.to_html(full_html=False, include_plotlyjs=True, default_width="760px")
    width, height = 1088, 624
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<style>
  body {{ font-family: Helvetica, Arial, sans-serif; margin: 0; width: {width}px;
         height: {height}px; color: #000000; background: #ffffff; }}
  .row {{ display: flex; align-items: flex-start; gap: 28px; padding: 12px 8px 0 0; }}
  aside {{ font-size: 17px; line-height: 1.4; padding-top: 28px; }}
  aside span {{ font-weight: 600; }}
</style>
</head>
<body>
<div class="row">
{body}
<aside>
{legend}
</aside>
</div>
</body>
</html>
"""
    html_path = PNG_OUT.with_suffix(".html")
    html_path.write_text(html)
    chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    subprocess.run(
        [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=2",
            f"--window-size={width},{height}",
            f"--screenshot={PNG_OUT}",
            "--virtual-time-budget=4000",
            html_path.as_uri(),
        ],
        check=True,
        capture_output=True,
    )
    html_path.unlink(missing_ok=True)
    print(f"wrote {PNG_OUT}")


def main() -> None:
    nucleosides = nucleoside_points()
    fig = scatter_figure(nucleosides, width=760, height=680)
    names = "".join(
        f"<div><span>{short}</span>: {name}</div>"
        for short, name in abbreviations(nucleosides)
    )
    body = fig.to_html(full_html=False, include_plotlyjs="cdn", default_width="760px")
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Reverse transcriptase Inhibitor Potency for HIV &amp; LINE-1</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Space+Grotesk:wght@500;600&display=swap" rel="stylesheet">
<style>
  body {{ font-family: Inter, sans-serif; margin: 24px auto; max-width: 1180px; color: #000000; background: #ffffff; }}
  h1 {{ font-family: "Space Grotesk", sans-serif; font-weight: 500; color: #000000; }}
  .row {{ display: flex; align-items: flex-start; gap: 36px; }}
  aside {{ color: #000000; font-size: 15px; line-height: 1.7; padding-top: 72px; }}
  aside span {{ font-weight: 600; }}
</style>
</head>
<body>
<h1>Reverse transcriptase Inhibitor Potency for HIV &amp; LINE-1</h1>
<div class="row">
{body}
<aside>
{names}
</aside>
</div>
<script>
document.fonts.ready.then(function () {{
  document.querySelectorAll(".js-plotly-plot").forEach(function (el) {{
    Plotly.relayout(el, {{}});
  }});
}});
</script>
</body>
</html>
"""
    OUT.write_text(html)
    fitted = [point for point in nucleosides if point.get("line1_censor") != "gt"]
    hiv = REFERENCE_NM / np.array([point["hiv_nM"] for point in fitted])
    line1 = REFERENCE_NM / np.array([point["line1_nM"] for point in fitted])
    fit_n = log_fit(hiv, line1)
    write_png(nucleosides)
    write_per_kb()
    print(f"wrote {OUT}")
    print(
        f"nucleosides n={fit_n['n']} slope={fit_n['slope']:.3f} "
        f"rho={fit_n['rho']:.3f} p={fit_n['p_rho']:.4g}"
    )
    for point in sorted(nucleosides, key=lambda item: item["short"].lower()):
        print(f"  {point['short']}: {point['name']}")


if __name__ == "__main__":
    main()
