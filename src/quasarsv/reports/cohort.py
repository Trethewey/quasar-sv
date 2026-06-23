"""Cohort dashboard — recurrence across samples, cohort circos overlay, per-sample matrix."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from ..model import FusionCall
from ..plots import circos_figure
from ..qc import cohort_recurrence
from ..classify import SOMATIC_CLINICAL, SOMATIC_WORTH_REVIEW, PHYSIOLOGICAL
from ..translocations import CANONICAL_CLASSES, group_by_translocation
from .common import (
    html_shell, fig_to_div, kpi_grid, render_table, tier_badge, Safe,
    gene_pair_with_badges,
)


def _per_sample_summary(calls: list[FusionCall], meta_index: dict | None = None) -> list[dict]:
    by_sample: dict[str, dict] = {}
    for c in calls:
        d = by_sample.setdefault(c.sample, {
            "tc": Counter(), "somatic": 0, "novel": 0, "physio": 0,
            "canonical": []})
        d["tc"][c.tier] += 1
        if c.event_class in SOMATIC_CLINICAL and c.tier in ("T1", "T2"):
            d["somatic"] += 1
            if c.event_class == "IG_driver_canonical":
                d["canonical"].append(f"{c.gene_a}-{c.gene_b}|{c.tier}")
        if c.event_class in SOMATIC_WORTH_REVIEW and c.tier in ("T1", "T2"):
            d["novel"] += 1
        if c.event_class in PHYSIOLOGICAL:
            d["physio"] += 1
    rows = []
    for s, d in sorted(by_sample.items()):
        tc = d["tc"]
        meta = (meta_index or {}).get(s)
        row = {
            "sample": Safe(f"<b>{s}</b>"),
            "cell_line": meta.cell_line if meta else "",
            "cohort": meta.cohort if meta else "",
            "ashm": meta.ashm_expected if meta else "",
            "cov": f"{meta.coverage:.0f}x" if meta and meta.coverage else "",
            "somatic": d["somatic"],
            "novel": d["novel"],
            "physio": d["physio"],
            "t1": tc.get("T1", 0),
            "t2": tc.get("T2", 0),
            "canonical_partners": Safe("; ".join(d["canonical"])) if d["canonical"] else "",
        }
        rows.append(row)
    return rows


def _recurrence_rows(calls: list[FusionCall], top: int = 50) -> list[dict]:
    rec = cohort_recurrence(calls, window=1000)
    if not rec:
        return []
    # Build a lookup from key -> sample example (gene names)
    by_key: dict[tuple, list[FusionCall]] = {}
    for c in calls:
        k = (c.chrom_a, c.pos_a // 1000, c.chrom_b, c.pos_b // 1000)
        by_key.setdefault(k, []).append(c)
    rows = []
    for k, n in rec.most_common(top):
        examples = by_key.get(k, [])
        if not examples:
            continue
        ex = max(examples, key=lambda c: c.n_callers)
        best_tier = next((t for t in ("T1", "T2", "T3") if any(c.tier == t for c in examples)), "T3")
        rows.append({
            "n_samples": n,
            "locus_a": Safe(f"<span class='mono'>{ex.gene_a or ex.chrom_a}:{ex.pos_a:,}</span>"),
            "locus_b": Safe(f"<span class='mono'>{ex.gene_b or ex.chrom_b}:{ex.pos_b:,}</span>"),
            "driver": ex.driver_locus,
            "known": Safe("<span class='badge partner'>known</span>") if ex.known_partner else "",
            "max_tier": tier_badge(best_tier),
        })
    return rows


def _best_per_sample_for_class(class_calls: list[FusionCall]) -> list[FusionCall]:
    """Return one best (lowest-tier, highest-evidence) call per sample for a class."""
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    best: dict[str, FusionCall] = {}
    for c in class_calls:
        cur = best.get(c.sample)
        score = (tier_rank.get(c.tier, 3), -(c.split_reads + c.discordant_pairs))
        if cur is None or score < (tier_rank.get(cur.tier, 3),
                                   -(cur.split_reads + cur.discordant_pairs)):
            best[c.sample] = c
    return list(best.values())


def _per_translocation_cards(calls: list[FusionCall]) -> str:
    grouped = group_by_translocation(calls)
    cards: list[str] = []
    for cls in CANONICAL_CLASSES:
        members = grouped.get(cls.key, [])
        if not members:
            continue
        best = _best_per_sample_for_class(members)
        best.sort(key=lambda c: ({"T1": 0, "T2": 1, "T3": 2}.get(c.tier, 3),
                                 -(c.split_reads + c.discordant_pairs)))
        rows = [{
            "sample": Safe(f"<b>{c.sample}</b>"),
            "tier_badge": tier_badge(c.tier),
            "gene_pair": gene_pair_with_badges(c.gene_a, c.gene_b),
            "coords": Safe(
                f"<span class='mono'>{c.chrom_a}:{c.pos_a:,}</span> &nbsp; "
                f"<span class='mono'>{c.chrom_b}:{c.pos_b:,}</span>"),
            "sr": c.split_reads, "pe": c.discordant_pairs,
            "callers": ", ".join(c.callers_supporting),
        } for c in best]
        cards.append(
            f"<div class='translocation-card'>"
            f"<h3>{cls.label} <span class='cyto'>{cls.cytoband}</span></h3>"
            f"<div class='disease'>{cls.disease} · "
            f"<b>{len(best)} sample{'s' if len(best) != 1 else ''}</b></div>"
            + render_table(rows, [
                ("Sample", "sample"), ("Tier", "tier_badge"),
                ("Genes", "gene_pair"), ("Coords", "coords"),
                ("SR", "sr"), ("PE", "pe"), ("Callers", "callers"),
            ])
            + "</div>"
        )
    if not cards:
        return "<p class='empty-row'>No canonical lymphoma translocations detected across the cohort.</p>"
    return "".join(cards)


def write_cohort_dashboard(
    calls: list[FusionCall],
    output_path: str,
    title: str = "Lymphoma fusion cohort dashboard",
    metadata_xlsx: str | None = None,
) -> str:
    nav_items = [
        ("Summary", "summary"),
        ("Circos", "circos"),
        ("Canonical translocations", "canon"),
        ("Per-sample", "samples"),
        ("Recurrence", "recurrence"),
    ]
    samples = sorted({c.sample for c in calls})
    tc = Counter(c.tier for c in calls)
    kpis = [
        ("samples", str(len(samples))),
        ("T1 total", str(tc.get("T1", 0))),
        ("T2 total", str(tc.get("T2", 0))),
        ("known partner T1+T2",
         str(sum(1 for c in calls if c.known_partner and c.tier in ("T1", "T2")))),
        ("recurrent (≥2 samples)", str(len(cohort_recurrence(calls)))),
    ]

    body = []
    body.append(f"<section class='section' id='summary'><h2>Summary</h2>{kpi_grid(kpis)}</section>")

    # Circos right under summary
    body.append("<section class='section' id='circos'>"
                "<h2>Cohort genome-wide view</h2>"
                "<p class='small'>Red arcs are canonical lymphoma translocations "
                "(labelled). Coloured ring segments mark the IG (blue) and TR "
                "(purple) loci, with intensity scaled by per-locus V(D)J / "
                "class-switch density.</p>"
                "<div class='plot-wrap'>"
                + fig_to_div(circos_figure(calls, title=f"cohort n={len(samples)}"))
                + "</div></section>")

    body.append("<section class='section' id='canon'>"
                "<h2>Canonical lymphoma translocations across the cohort</h2>"
                "<p class='small'>One card per recurrent translocation class. "
                "Each card lists every sample where that translocation was "
                "detected (best evidence per sample). Physiological IG V(D)J / "
                "class-switch signal is excluded.</p>"
                + _per_translocation_cards(calls)
                + "</section>")

    meta_index = {}
    if metadata_xlsx:
        try:
            from ..metadata import load_cohort_metadata_xlsx, build_metadata_index
            meta_index = build_metadata_index(load_cohort_metadata_xlsx(metadata_xlsx))
        except Exception:
            meta_index = {}

    cols = [("Sample", "sample"), ("Cell line", "cell_line"),
            ("Cohort", "cohort"), ("aSHM", "ashm"), ("Cov.", "cov"),
            ("Somatic T1+T2", "somatic"),
            ("Novel T1+T2", "novel"),
            ("Physio. IG/TR", "physio"),
            ("Tier T1", "t1"), ("Tier T2", "t2"),
            ("Canonical partners", "canonical_partners")]
    body.append("<section class='section' id='samples'>"
                "<h2>Per-sample summary</h2>"
                "<p class='small'>Somatic = IG-driver canonical + driver-driver "
                "(clinically actionable). Novel = IG-driver non-canonical or "
                "driver-intergenic (worth orthogonal review). "
                "Physio. = intra-IG/TR V(D)J + class-switch (not pathological).</p>"
                + render_table(_per_sample_summary(calls, meta_index), cols)
                + "</section>")

    body.append("<section class='section' id='recurrence'><h2>Recurrent rearrangements</h2>"
                "<p class='small'>Breakpoint pairs seen in ≥2 samples (±1 kb).</p>"
                + render_table(_recurrence_rows(calls), [
                    ("n samples", "n_samples"),
                    ("Locus A", "locus_a"), ("Locus B", "locus_b"),
                    ("Driver", "driver"), ("Known", "known"),
                ]) + "</section>")

    html = html_shell(title, nav_items, "\n".join(body))
    Path(output_path).write_text(html, encoding="utf-8")
    return output_path
