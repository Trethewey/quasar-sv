"""Locus-level gene-rearrangement plots — show breakpoints landing in or near
specified driver loci with mate annotations.
"""
from __future__ import annotations

from collections import defaultdict

import plotly.graph_objects as go

from ..model import FusionCall
from ..annotate import GeneEntry, load_builtin_loci


def _matches_locus(c: FusionCall, locus: GeneEntry, pad: int = 50_000) -> str | None:
    """Return 'a' or 'b' if a breakpoint of `c` falls in the locus (with pad)."""
    if c.chrom_a == locus.chrom and locus.start - pad <= c.pos_a <= locus.end + pad:
        return "a"
    if c.chrom_b == locus.chrom and locus.start - pad <= c.pos_b <= locus.end + pad:
        return "b"
    return None


def locus_figure(
    calls: list[FusionCall],
    gene: str,
    pad: int = 50_000,
) -> go.Figure:
    """Track plot showing all breakpoints near a locus with partner labels."""
    loci = {g.gene: g for g in load_builtin_loci()}
    locus = loci.get(gene)
    if locus is None:
        return _empty(f"Unknown locus {gene}")
    hits: list[tuple[FusionCall, str]] = []
    for c in calls:
        side = _matches_locus(c, locus, pad)
        if side:
            hits.append((c, side))
    if not hits:
        return _empty(f"No calls near {gene} (chr{locus.chrom}:{locus.start}-{locus.end})")

    fig = go.Figure()
    # Gene body
    fig.add_shape(type="rect",
                  x0=locus.start, x1=locus.end, y0=-0.2, y1=0.2,
                  fillcolor="#16213e", line=dict(color="#16213e"),
                  layer="below")
    fig.add_annotation(x=(locus.start + locus.end) / 2, y=0.45,
                       text=f"<b>{locus.gene}</b> ({locus.role})",
                       showarrow=False, font=dict(size=12, color="#16213e"))

    # Breakpoint markers, lane = tier
    lane = {"T1": 1.6, "T2": 1.0, "T3": 0.5}
    color = {"T1": "#d62728", "T2": "#ff7f0e", "T3": "#9ecae1"}
    for c, side in hits:
        x = c.pos_a if side == "a" else c.pos_b
        y = lane.get(c.tier, 0.5)
        partner = (
            f"{c.gene_b or c.chrom_b}:{c.pos_b}" if side == "a"
            else f"{c.gene_a or c.chrom_a}:{c.pos_a}"
        )
        col = color.get(c.tier, "#888")
        fig.add_trace(go.Scatter(
            x=[x], y=[y], mode="markers",
            marker=dict(size=max(6, min(18, 6 + c.split_reads * 0.4)),
                        color=col, line=dict(color="white", width=1)),
            hovertext=(
                f"<b>{partner}</b><br>{c.tier} | callers={','.join(c.callers_supporting)}"
                f"<br>SR={c.split_reads} PE={c.discordant_pairs} AS={c.assembly_contigs}"
                + (f"<br>known: {c.known_partner_source}" if c.known_partner else "")
            ),
            hoverinfo="text", showlegend=False,
        ))

    # X axis padded
    x_min = min(locus.start, min(_x(c, s) for c, s in hits)) - 5000
    x_max = max(locus.end, max(_x(c, s) for c, s in hits)) + 5000
    fig.update_layout(
        title=f"Rearrangements at {gene} (chr{locus.chrom}:{locus.start:,}-{locus.end:,})",
        xaxis=dict(title=f"chr{locus.chrom} (bp)", range=[x_min, x_max], showgrid=False),
        yaxis=dict(showgrid=False, range=[-0.7, 2.0],
                   tickvals=[0.5, 1.0, 1.6], ticktext=["T3", "T2", "T1"]),
        height=300, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=60, r=20, t=60, b=40),
    )
    return fig


def _x(c: FusionCall, side: str) -> int:
    return c.pos_a if side == "a" else c.pos_b


def locus_summary_table(
    calls: list[FusionCall],
    genes: list[str] | None = None,
    pad: int = 50_000,
) -> list[dict]:
    """Per-locus call count summary suitable for HTML tabulation."""
    loci = load_builtin_loci()
    if genes:
        loci = [g for g in loci if g.gene in set(genes)]
    by_locus: dict[str, dict] = {}
    for g in loci:
        by_locus[g.gene] = {"gene": g.gene, "chrom": g.chrom,
                            "start": g.start, "end": g.end, "role": g.role,
                            "T1": 0, "T2": 0, "T3": 0, "known_partner": 0}
    for c in calls:
        for g in loci:
            if _matches_locus(c, g, pad):
                by_locus[g.gene][c.tier] = by_locus[g.gene].get(c.tier, 0) + 1
                if c.known_partner:
                    by_locus[g.gene]["known_partner"] += 1
    return [v for v in by_locus.values() if v["T1"] + v["T2"] + v["T3"] > 0]


def _empty(msg: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
    fig.update_layout(height=200, paper_bgcolor="white", plot_bgcolor="white",
                      xaxis=dict(visible=False), yaxis=dict(visible=False),
                      margin=dict(l=20, r=20, t=20, b=20))
    return fig
