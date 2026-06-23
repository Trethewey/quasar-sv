"""Validation plots — per-caller concordance and (when truth available) ROC.

Without a labelled truth set we report:
  * concordance heatmap: pairwise Jaccard between callers
  * evidence breakdown: how much each evidence type contributes per caller
  * replicate concordance: scatter of T1+T2 calls in replicate pairs
"""
from __future__ import annotations

from collections import Counter, defaultdict

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..model import FusionCall


def _bp_key(c: FusionCall, window: int) -> tuple:
    return (c.chrom_a, c.pos_a // window, c.chrom_b, c.pos_b // window,
            c.strand_a, c.strand_b)


def per_caller_concordance(calls: list[FusionCall], window: int = 500) -> go.Figure:
    """Pairwise Jaccard concordance heatmap across callers."""
    callers_used: dict[str, set] = defaultdict(set)
    for c in calls:
        k = _bp_key(c, window)
        for cl in c.callers_supporting:
            callers_used[cl].add(k)
    callers = sorted(callers_used.keys())
    n = len(callers)
    if n == 0:
        return _empty_fig("No caller-tagged calls to compare")
    mat = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(callers):
        for j, b in enumerate(callers):
            sa = callers_used[a]
            sb = callers_used[b]
            inter = len(sa & sb)
            uni = len(sa | sb)
            mat[i][j] = inter / uni if uni else 0.0
    fig = go.Figure(data=go.Heatmap(
        z=mat, x=callers, y=callers,
        colorscale="Blues", zmin=0, zmax=1,
        text=[[f"{v:.2f}" for v in row] for row in mat],
        texttemplate="%{text}",
        hovertemplate="<b>%{y}</b> vs <b>%{x}</b><br>Jaccard %{z:.3f}<extra></extra>",
        colorbar=dict(title="Jaccard"),
    ))
    fig.update_layout(
        title="Per-caller breakpoint concordance (Jaccard)",
        height=380,
        margin=dict(l=60, r=20, t=60, b=40),
        paper_bgcolor="white", plot_bgcolor="white",
    )
    return fig


def per_caller_evidence_breakdown(calls: list[FusionCall]) -> go.Figure:
    """Stacked bar: per-caller fraction of T1 / T2 / T3 calls."""
    rows: dict[str, Counter] = defaultdict(Counter)
    for c in calls:
        for cl in c.callers_supporting:
            rows[cl][c.tier] += 1
    callers = sorted(rows.keys())
    if not callers:
        return _empty_fig("No callers present")
    fig = go.Figure()
    for tier, colour in (("T1", "#d62728"), ("T2", "#ff7f0e"), ("T3", "#9ecae1")):
        fig.add_trace(go.Bar(
            x=callers, y=[rows[cl].get(tier, 0) for cl in callers],
            name=tier, marker_color=colour,
        ))
    fig.update_layout(
        barmode="stack",
        title="Per-caller tier contribution",
        height=350,
        margin=dict(l=60, r=20, t=60, b=40),
        paper_bgcolor="white", plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
    )
    fig.update_yaxes(title="Number of fusion candidates")
    return fig


def replicate_concordance(
    calls_by_sample: dict[str, list[FusionCall]],
    replicate_pairs: list[tuple[str, str]],
    window: int = 500,
) -> go.Figure:
    """Scatter: shared T1/T2 calls between replicate sample pairs.

    `replicate_pairs` is a list of (sample_a, sample_b). Plots per-pair shared
    versus pair-unique counts (T1/T2 only).
    """
    if not replicate_pairs:
        return _empty_fig("No replicate pairs supplied")
    points = []
    for a, b in replicate_pairs:
        ca = {_bp_key(c, window) for c in calls_by_sample.get(a, []) if c.tier in ("T1", "T2")}
        cb = {_bp_key(c, window) for c in calls_by_sample.get(b, []) if c.tier in ("T1", "T2")}
        shared = len(ca & cb)
        a_only = len(ca - cb)
        b_only = len(cb - ca)
        points.append((a, b, shared, a_only, b_only))

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Shared T1/T2 calls per pair",
                                        "Jaccard per pair"))
    labels = [f"{p[0][:12]}—{p[1][:12]}" for p in points]
    shared = [p[2] for p in points]
    jaccard = [p[2] / (p[2] + p[3] + p[4]) if (p[2] + p[3] + p[4]) else 0 for p in points]
    fig.add_trace(go.Bar(x=labels, y=shared, marker_color="#16213e"), row=1, col=1)
    fig.add_trace(go.Bar(x=labels, y=jaccard, marker_color="#d62728"), row=1, col=2)
    fig.update_layout(
        height=380, showlegend=False, paper_bgcolor="white",
        margin=dict(l=60, r=20, t=60, b=80),
    )
    fig.update_xaxes(tickangle=-30)
    return fig


def precision_recall(calls: list[FusionCall], truth_keys: set[tuple]) -> go.Figure:
    """ROC-like PR curve as a function of split-read threshold.

    `truth_keys` is a set of (chrom_a, pos_a//window, chrom_b, pos_b//window, ...).
    Only used when a labelled truth set is provided.
    """
    if not truth_keys:
        return _empty_fig("No truth set supplied — precision/recall unavailable")
    thresholds = list(range(0, 21))
    prec = []
    rec = []
    for thr in thresholds:
        positives = [c for c in calls if c.split_reads >= thr]
        if not positives:
            prec.append(1.0)
            rec.append(0.0)
            continue
        tp = sum(1 for c in positives if _bp_key(c, 500) in truth_keys)
        prec.append(tp / len(positives))
        rec.append(tp / max(len(truth_keys), 1))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=rec, y=prec, mode="lines+markers",
                             line=dict(color="#16213e", width=2)))
    fig.update_layout(
        title="Precision / recall vs split-read threshold",
        xaxis=dict(title="Recall", range=[0, 1]),
        yaxis=dict(title="Precision", range=[0, 1.02]),
        height=360, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=60, r=20, t=60, b=40),
    )
    return fig


def _empty_fig(msg: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
    fig.update_layout(height=200, paper_bgcolor="white", plot_bgcolor="white",
                      xaxis=dict(visible=False), yaxis=dict(visible=False),
                      margin=dict(l=20, r=20, t=20, b=20))
    return fig
