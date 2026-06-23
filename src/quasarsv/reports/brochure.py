"""Per-sample clinical brochure — tables-first, then locus close-ups, then circos."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from ..model import FusionCall
from ..plots import (
    circos_figure, qc_summary_figure, evidence_distribution,
    locus_figure, locus_summary_table,
)
from ..classify import (
    SOMATIC_CLINICAL, SOMATIC_WORTH_REVIEW, PHYSIOLOGICAL, is_somatic_clinical,
)
from ..translocations import (
    CANONICAL_CLASSES, classify_translocation, group_by_translocation, is_ig_involved,
)
from .common import (
    html_shell, fig_to_div, kpi_grid, render_table, tier_badge, Safe,
    gene_pair_with_badges, locus_badge,
)


CLINICAL_LOCI = ["MYC", "BCL2", "BCL6", "CCND1", "IRF4", "DUSP22", "ALK", "MALT1", "BCL10",
                 "IGH", "IGK", "IGL", "PAX5", "REL", "BCL3"]


def _is_rescue(call: FusionCall) -> bool:
    return "inferred_via_artefact_rescue" in call.qc_flags


def _kpis(calls: list[FusionCall]) -> list[tuple[str, str]]:
    tc = Counter(c.tier for c in calls)
    n_somatic = sum(1 for c in calls if c.event_class in SOMATIC_CLINICAL
                    and c.tier in ("T1", "T2"))
    n_novel = sum(1 for c in calls if c.event_class in SOMATIC_WORTH_REVIEW
                  and c.tier in ("T1", "T2"))
    n_physio = sum(1 for c in calls if c.event_class in PHYSIOLOGICAL)
    n_callers = sorted({cl for c in calls for cl in c.callers_supporting})
    return [
        ("Somatic (T1+T2)", str(n_somatic)),
        ("Novel/review (T1+T2)", str(n_novel)),
        ("Physiological IG/TR", str(n_physio)),
        ("T1 total", str(tc.get("T1", 0))),
        ("T2 total", str(tc.get("T2", 0))),
        ("Callers run", str(len(n_callers))),
    ]


def _t1_t2_table(calls: list[FusionCall], include_rescue: bool = False) -> list[dict]:
    rows = []
    for c in sorted(calls, key=lambda x: (
            {"T1": 0, "T2": 1, "T3": 2}.get(x.tier, 3),
            -x.split_reads, -x.n_callers)):
        if c.tier == "T3":
            continue
        if _is_rescue(c) != include_rescue:
            continue
        rows.append({
            "tier_badge": tier_badge(c.tier),
            "locus_a": Safe(f"<span class='mono'>{c.gene_a or c.chrom_a}:{c.pos_a:,}{c.strand_a}</span>"),
            "locus_b": Safe(f"<span class='mono'>{c.gene_b or c.chrom_b}:{c.pos_b:,}{c.strand_b}</span>"),
            "sv_type": c.sv_type,
            "callers": ", ".join(c.callers_supporting),
            "n_ev": c.n_evidence_types,
            "sr": c.split_reads,
            "pe": c.discordant_pairs,
            "as": c.assembly_contigs,
            "vaf": f"{c.vaf:.3f}" if c.vaf else "",
            "driver": c.driver_locus or "",
            "known": Safe("<span class='badge partner'>known</span>") if c.known_partner else "",
            "flags": Safe(" ".join(f"<span class='badge warn'>{f}</span>" for f in c.qc_flags)),
        })
    return rows


def _driver_locus_table(calls: list[FusionCall]) -> list[dict]:
    rows = locus_summary_table(calls, genes=CLINICAL_LOCI)
    rows.sort(key=lambda r: (-r["T1"], -r["T2"]))
    return [{
        "gene": Safe(f"<b>{r['gene']}</b>"),
        "role": r["role"],
        "coords": Safe(f"<span class='mono'>chr{r['chrom']}:{r['start']:,}-{r['end']:,}</span>"),
        "t1": r["T1"], "t2": r["T2"], "t3": r["T3"],
        "known_partner": r["known_partner"],
    } for r in rows]


def write_brochure(
    sample: str,
    calls: list[FusionCall],
    output_path: str,
    title: str | None = None,
    metadata=None,
) -> str:
    title = title or f"Fusion detection report — {sample}"
    nav_items = [
        ("Summary", "summary"),
        ("Circos", "circos"),
        ("Canonical somatic", "somatic"),
        ("Driver ↔ driver", "driver-driver"),
        ("Novel (review)", "novel"),
        ("Physiological IG/TR", "physio"),
        ("Driver-locus hits", "drivers"),
        ("Locus close-ups", "loci"),
        ("QC", "qc"),
    ]

    body_parts: list[str] = []
    meta_html = ""
    if metadata:
        meta_kpis = [
            ("cell line", metadata.cell_line or "—"),
            ("cohort", metadata.cohort or "—"),
            ("aSHM expected", metadata.ashm_expected or "—"),
            ("coverage", f"{metadata.coverage:.0f}×" if metadata.coverage else "—"),
            ("project", metadata.project or "—"),
            ("run", metadata.run or "—"),
        ]
        meta_html = (
            "<section class='section'><h2>Sample provenance</h2>"
            + kpi_grid(meta_kpis)
            + "</section>"
        )
    body_parts.append(meta_html)
    body_parts.append(
        "<section class='section' id='summary'><h2>Detection summary</h2>"
        + kpi_grid(_kpis(calls))
        + "</section>"
    )

    # Move circos to right under summary
    body_parts.append(
        "<section class='section' id='circos'>"
        "<h2>Genome-wide view</h2>"
        "<p class='small'>Red arcs are canonical lymphoma translocations, "
        "labelled with the partner gene names. Pale blue arcs in the "
        "background are physiological IG / TR V(D)J + class-switch "
        "rearrangements.</p>"
        "<div class='plot-wrap'>"
        + fig_to_div(circos_figure(calls, title=sample))
        + "</div></section>"
    )

    tier_rank = {"T1": 0, "T2": 1, "T3": 2}

    def _row(c):
        return {
            "tier_badge": tier_badge(c.tier),
            "gene_pair": gene_pair_with_badges(c.gene_a, c.gene_b),
            "coords": Safe(
                f"<span class='mono'>{c.chrom_a}:{c.pos_a:,}</span> &nbsp; "
                f"<span class='mono'>{c.chrom_b}:{c.pos_b:,}</span>"),
            "sr": c.split_reads, "pe": c.discordant_pairs,
            "callers": ", ".join(c.callers_supporting),
            "qc": Safe(" ".join(f"<span class='badge warn'>{f}</span>" for f in c.qc_flags)),
        }

    # ── Per-canonical-translocation cards ──
    grouped = group_by_translocation(calls)
    canon_cards: list[str] = []
    n_canon_total = 0
    for cls in CANONICAL_CLASSES:
        members = grouped.get(cls.key, [])
        if not members:
            continue
        members.sort(key=lambda x: (tier_rank.get(x.tier, 3),
                                    -(x.split_reads + x.discordant_pairs)))
        rows = [_row(c) for c in members]
        n_canon_total += len(members)
        canon_cards.append(
            f"<div class='translocation-card'>"
            f"<h3>{cls.label} <span class='cyto'>{cls.cytoband}</span></h3>"
            f"<div class='disease'>{cls.disease}</div>"
            + render_table(rows, [
                ("Tier", "tier_badge"),
                ("Genes", "gene_pair"),
                ("Coords", "coords"),
                ("SR", "sr"), ("PE", "pe"),
                ("Callers", "callers"),
                ("Flags", "qc"),
            ])
            + "</div>"
        )
    body_parts.append(
        "<section class='section' id='somatic'>"
        "<h2>Canonical somatic translocations</h2>"
        "<p class='small'>One card per recurrent lymphoma translocation class. "
        "Calls only appear here if at least one breakpoint matches the canonical "
        f"gene pair. Total canonical calls in this sample: <b>{n_canon_total}</b>.</p>"
        + ("".join(canon_cards) if canon_cards
           else "<p class='empty-row'>No canonical lymphoma translocations detected in this sample.</p>")
        + "</section>"
    )

    # ── Driver-driver fusions (non-IG, non-canonical-class) ──
    drv_drv_rows = [
        _row(c) for c in sorted(
            (c for c in calls
             if c.event_class == "driver_driver"
             and classify_translocation(c) is None
             and c.tier in ("T1", "T2")),
            key=lambda x: (tier_rank.get(x.tier, 3), -(x.split_reads + x.discordant_pairs)),
        )
    ]
    body_parts.append(
        "<section class='section' id='driver-driver'>"
        "<h2>Driver ↔ driver fusions (non-canonical, T1+T2)</h2>"
        "<p class='small'>Two annotated driver genes joined together (not "
        "matching a known recurrent translocation class).</p>"
        + render_table(drv_drv_rows, [
            ("Tier", "tier_badge"), ("Genes", "gene_pair"),
            ("Coords", "coords"), ("SR", "sr"), ("PE", "pe"),
            ("Callers", "callers"), ("Flags", "qc"),
        ])
        + "</section>"
    )

    # ── Novel / review-worthy ──
    novel_rows = [
        _row(c) for c in sorted(
            (c for c in calls if c.event_class in SOMATIC_WORTH_REVIEW
             and c.tier in ("T1", "T2")),
            key=lambda x: (tier_rank.get(x.tier, 3), -(x.split_reads + x.discordant_pairs)),
        )[:60]
    ]
    body_parts.append(
        "<section class='section' id='novel'>"
        "<h2>Putative novel / non-canonical (T1+T2)</h2>"
        "<p class='small'>Driver rearranged to an unannotated region, or IG ↔ "
        "driver pair not in any canonical class. Orthogonal review recommended. "
        "Top 60 by evidence shown.</p>"
        + render_table(novel_rows, [
            ("Tier", "tier_badge"), ("Genes", "gene_pair"),
            ("Coords", "coords"), ("SR", "sr"), ("PE", "pe"),
            ("Callers", "callers"), ("Flags", "qc"),
        ])
        + "</section>"
    )

    # ── Physiological IG/TR (summarised, not enumerated) ──
    physio = [c for c in calls if c.event_class in PHYSIOLOGICAL]
    physio_counter = Counter(
        f"{c.gene_a} (intra)" if c.gene_a == c.gene_b
        else f"{c.gene_a} ↔ {c.gene_b}"
        for c in physio
    )
    physio_summary_rows = [
        {"locus": Safe(locus_badge(k.split(" ")[0])
                       + (Safe(" intra")) if "(intra)" in k
                       else gene_pair_with_badges(*k.split(" ↔ "))),
         "count": v} for k, v in physio_counter.most_common(15)
    ]
    body_parts.append(
        "<section class='section' id='physio'>"
        "<h2>Physiological IG / TR rearrangements</h2>"
        "<p class='small'>Intra-IG and intra-TR breakpoints from V(D)J "
        "recombination and class-switch recombination — expected in any "
        f"B- or T-cell sample. <b>{len(physio)}</b> total physiological calls in "
        "this sample (top 15 loci shown).</p>"
        + render_table(physio_summary_rows, [
            ("Locus", "locus"), ("Count", "count"),
        ])
        + "</section>"
    )

    body_parts.append(
        "<section class='section' id='drivers'>"
        "<h2>Driver-locus hits</h2>"
        "<p class='small'>Breakpoint counts within ±50 kb of each lymphoma driver locus.</p>"
        + render_table(_driver_locus_table(calls), [
            ("Gene", "gene"), ("Role", "role"), ("Coords", "coords"),
            ("T1", "t1"), ("T2", "t2"), ("T3", "t3"),
            ("Known partner", "known_partner"),
        ])
        + "</section>"
    )


    loci_html: list[str] = []
    for g in CLINICAL_LOCI:
        hits = [c for c in calls
                if (c.gene_a == g or c.gene_b == g or c.driver_locus == g)]
        if not hits:
            continue
        loci_html.append(f"<h3>{g}</h3>")
        loci_html.append("<div class='plot-wrap'>" + fig_to_div(locus_figure(calls, g)) + "</div>")
    body_parts.append(
        "<section class='section' id='loci'><h2>Locus close-ups</h2>"
        + ("".join(loci_html) if loci_html else "<p class='small'>No clinical driver loci hit.</p>")
        + "</section>"
    )

    body_parts.append(
        "<section class='section' id='qc'><h2>QC summary</h2>"
        "<div class='plot-wrap'>"
        + fig_to_div(qc_summary_figure(calls, sample=sample))
        + "</div>"
        "<div class='plot-wrap'>"
        + fig_to_div(evidence_distribution(calls))
        + "</div></section>"
    )

    html = html_shell(title, nav_items, "\n".join(body_parts))
    Path(output_path).write_text(html, encoding="utf-8")
    return output_path
