"""Validation report — per-caller performance, concordance, replicate stability."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from ..model import FusionCall
from ..plots.validation import (
    per_caller_concordance, per_caller_evidence_breakdown,
    replicate_concordance, precision_recall,
)
from ..plots.qc_plots import qc_summary_figure, evidence_distribution
from .common import html_shell, fig_to_div, kpi_grid


def write_validation_report(
    calls: list[FusionCall],
    output_path: str,
    replicate_pairs: list[tuple[str, str]] | None = None,
    truth_keys: set[tuple] | None = None,
    title: str = "quasarsv validation report",
) -> str:
    nav_items = [
        ("Caller agreement", "concordance"),
        ("Tier composition", "composition"),
        ("Replicates", "reps"),
        ("Precision/recall", "pr"),
        ("Evidence", "ev"),
    ]
    callers = sorted({cl for c in calls for cl in c.callers_supporting})
    kpis = [
        ("callers", str(len(callers))),
        ("samples", str(len({c.sample for c in calls}))),
        ("T1 calls", str(sum(1 for c in calls if c.tier == "T1"))),
        ("T2 calls", str(sum(1 for c in calls if c.tier == "T2"))),
        ("replicate pairs", str(len(replicate_pairs or []))),
        ("truth set", "supplied" if truth_keys else "—"),
    ]

    by_sample: dict[str, list[FusionCall]] = defaultdict(list)
    for c in calls:
        by_sample[c.sample].append(c)

    body = []
    body.append(f"<section class='section'><h2>Header</h2>{kpi_grid(kpis)}</section>")

    body.append("<section class='section' id='concordance'>"
                "<h2>Per-caller agreement</h2>"
                "<p class='small'>Jaccard similarity of breakpoints (±500 bp).</p>"
                "<div class='plot-wrap'>"
                + fig_to_div(per_caller_concordance(calls))
                + "</div></section>")

    body.append("<section class='section' id='composition'>"
                "<h2>Tier composition per caller</h2>"
                "<div class='plot-wrap'>"
                + fig_to_div(per_caller_evidence_breakdown(calls))
                + "</div></section>")

    body.append("<section class='section' id='reps'>"
                "<h2>Replicate concordance</h2>"
                "<div class='plot-wrap'>"
                + fig_to_div(replicate_concordance(by_sample, replicate_pairs or []))
                + "</div></section>")

    body.append("<section class='section' id='pr'>"
                "<h2>Precision / recall vs split-read threshold</h2>"
                "<div class='plot-wrap'>"
                + fig_to_div(precision_recall(calls, truth_keys or set()))
                + "</div></section>")

    body.append("<section class='section' id='ev'>"
                "<h2>QC and evidence distributions</h2>"
                "<div class='plot-wrap'>"
                + fig_to_div(qc_summary_figure(calls))
                + "</div>"
                "<div class='plot-wrap'>"
                + fig_to_div(evidence_distribution(calls))
                + "</div></section>")

    html = html_shell(title, nav_items, "\n".join(body))
    Path(output_path).write_text(html, encoding="utf-8")
    return output_path
