"""QC plots — tier histograms, evidence-type distributions, sample-level metrics."""
from __future__ import annotations

from collections import Counter

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..model import FusionCall


def qc_summary_figure(calls: list[FusionCall], sample: str | None = None) -> go.Figure:
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("Tier distribution", "Caller support", "QC flags"),
        specs=[[{"type": "bar"}, {"type": "bar"}, {"type": "bar"}]],
    )
    tier_c = Counter(c.tier for c in calls)
    fig.add_trace(go.Bar(
        x=["T1", "T2", "T3"],
        y=[tier_c.get(t, 0) for t in ("T1", "T2", "T3")],
        marker_color=["#d62728", "#ff7f0e", "#9ecae1"],
        showlegend=False,
    ), row=1, col=1)

    caller_c = Counter()
    for c in calls:
        caller_c[len(c.callers_supporting)] += 1
    n_max = max(caller_c.keys(), default=0)
    fig.add_trace(go.Bar(
        x=[str(k) for k in range(1, n_max + 1)],
        y=[caller_c.get(k, 0) for k in range(1, n_max + 1)],
        marker_color="#16213e", showlegend=False,
    ), row=1, col=2)

    flag_c = Counter()
    for c in calls:
        for f in c.qc_flags:
            flag_c[f] += 1
    flags = list(flag_c.keys()) or ["(none)"]
    counts = [flag_c.get(f, 0) for f in flags]
    fig.add_trace(go.Bar(x=flags, y=counts, marker_color="#7f7f7f",
                        showlegend=False), row=1, col=3)

    title = f"QC summary — {sample}" if sample else "QC summary"
    fig.update_layout(
        title=title, height=360,
        paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=60, r=20, t=60, b=40),
    )
    return fig


def evidence_distribution(calls: list[FusionCall]) -> go.Figure:
    """Histogram of split-read / discordant-pair / assembly-contig support."""
    sr = [c.split_reads for c in calls if c.split_reads > 0]
    pe = [c.discordant_pairs for c in calls if c.discordant_pairs > 0]
    asm = [c.assembly_contigs for c in calls if c.assembly_contigs > 0]
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=sr, name="split reads",
                               marker_color="#d62728", opacity=0.7, nbinsx=40))
    fig.add_trace(go.Histogram(x=pe, name="discordant pairs",
                               marker_color="#16213e", opacity=0.7, nbinsx=40))
    fig.add_trace(go.Histogram(x=asm, name="assembly contigs",
                               marker_color="#ff7f0e", opacity=0.7, nbinsx=40))
    fig.update_layout(
        title="Evidence-type support per fusion candidate",
        barmode="overlay", height=320,
        xaxis=dict(title="evidence count"),
        yaxis=dict(title="number of candidates", type="log"),
        paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=60, r=20, t=60, b=40),
        legend=dict(orientation="h", y=1.05),
    )
    return fig
