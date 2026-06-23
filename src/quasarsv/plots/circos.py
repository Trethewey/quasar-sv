"""Circos / chord-style genome-wide rearrangement plot using Plotly polar.

Each chromosome is a sector; each rearrangement is an arc between
breakpoint angles. Colour by tier, width by split-read support.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import plotly.graph_objects as go

from ..model import FusionCall

# GRCh38 chromosome lengths (bp) — bare-numeric names; X/Y included.
CHROM_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}
ORDER = list(CHROM_LENGTHS.keys())

TIER_COLOR = {"T1": "#d62728", "T2": "#ff7f0e", "T3": "#9ecae1"}


@dataclass
class ChromSector:
    chrom: str
    angle_start: float
    angle_end: float
    length: int


def _build_sectors(gap_deg: float = 0.6) -> dict[str, ChromSector]:
    total = sum(CHROM_LENGTHS.values())
    free = 360 - gap_deg * len(ORDER)
    out: dict[str, ChromSector] = {}
    a = 0.0
    for c in ORDER:
        ln = CHROM_LENGTHS[c]
        span = (ln / total) * free
        out[c] = ChromSector(c, a, a + span, ln)
        a += span + gap_deg
    return out


def _pos_to_angle(chrom: str, pos: int, sectors: dict[str, ChromSector]) -> float | None:
    s = sectors.get(_strip_chr(chrom))
    if s is None:
        return None
    frac = min(max(pos / s.length, 0.0), 1.0)
    deg = s.angle_start + frac * (s.angle_end - s.angle_start)
    return math.radians(deg)


def _strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def _arc_points(theta1: float, theta2: float, r: float, n: int = 50):
    # quadratic Bezier from (r,theta1) to (r,theta2) bending through origin (chord interior)
    x1, y1 = r * math.cos(theta1), r * math.sin(theta1)
    x2, y2 = r * math.cos(theta2), r * math.sin(theta2)
    # control point pulled toward centre
    cx = (x1 + x2) * 0.15
    cy = (y1 + y2) * 0.15
    xs, ys = [], []
    for i in range(n + 1):
        t = i / n
        x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cx + t ** 2 * x2
        y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cy + t ** 2 * y2
        # convert back to polar
        rr = math.hypot(x, y)
        tt = math.atan2(y, x)
        xs.append(rr)
        ys.append(math.degrees(tt) % 360)
    return ys, xs   # (theta_deg, r)


def circos_figure(
    calls: list[FusionCall],
    title: str = "",
    height: int = 720,
) -> go.Figure:
    """Polar circos figure.

    Renders **only** two categories:
      * canonical lymphoma translocations (red, labelled with partner pair)
      * physiological IG/TR rearrangements (very pale background — gives a
        sense of V(D)J + class-switch density without dominating the plot)

    Non-canonical, non-IG calls are excluded — they belong in their own
    tables, not on the genome ring.
    """
    from ..translocations import classify_translocation, is_ig_involved
    from ..classify import PHYSIOLOGICAL

    sectors = _build_sectors()
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    fig = go.Figure()
    R_OUT = 1.0
    R_IN = 0.92

    # chromosome ring
    for c, s in sectors.items():
        thetas = list(_linspace(s.angle_start, s.angle_end, 60))
        r_out = [R_OUT] * len(thetas)
        fig.add_trace(go.Scatterpolar(
            theta=thetas + thetas[::-1],
            r=r_out + [R_IN] * len(thetas),
            mode="lines",
            fill="toself",
            line=dict(width=0),
            fillcolor="#16213e",
            hovertext=f"chr{c}",
            hoverinfo="text",
            showlegend=False,
        ))
        # Chromosome label sits ON the dark ring (white text)
        mid = (s.angle_start + s.angle_end) / 2
        fig.add_trace(go.Scatterpolar(
            theta=[mid], r=[0.96],
            mode="text", text=[f"<b>{c}</b>"],
            textfont=dict(size=13, color="white"),
            showlegend=False, hoverinfo="skip",
        ))

    # ── Physiological IG/TR — highlight the IG/TR loci on the ring ──
    # Intra-locus arcs collapse to invisible points, so instead we paint a
    # thick coloured ring segment at each IG/TR locus, with opacity scaled by
    # the per-locus breakpoint count. Inter-IG arcs (IGH↔IGK etc., rare but
    # real) are drawn as faint blue arcs.
    physio = [c for c in calls if c.event_class in PHYSIOLOGICAL]
    intra_counts: dict[str, int] = {}
    inter_pairs: list = []
    for c in physio:
        if c.event_class == "IG_intra":
            intra_counts[c.gene_a or c.gene_b] = intra_counts.get(c.gene_a or c.gene_b, 0) + 1
        else:
            inter_pairs.append(c)

    # IG/TR ring highlights
    from ..annotate import load_builtin_loci
    loci = {g.gene: g for g in load_builtin_loci()}
    IG_FILL = {  # darker = B-cell IG (cooler), purple-ish = TR (T-cell)
        "IGH": "#1d4ed8", "IGK": "#1d4ed8", "IGL": "#1d4ed8",
        "IGH_Emu": "#1d4ed8", "IGH_3RR": "#1d4ed8",
        "TRA": "#7c3aed", "TRB": "#7c3aed",
        "TRG": "#7c3aed", "TRD": "#7c3aed",
    }
    max_count = max(intra_counts.values(), default=1)
    R_HL_OUT = 1.18
    R_HL_IN = 1.08
    for gene, count in intra_counts.items():
        g = loci.get(gene)
        if g is None:
            continue
        a_start_rad = _pos_to_angle(g.chrom, g.start, sectors)
        a_end_rad = _pos_to_angle(g.chrom, g.end, sectors)
        if a_start_rad is None or a_end_rad is None:
            continue
        a_start = math.degrees(a_start_rad)
        a_end = math.degrees(a_end_rad)
        # Make the IG/TR band clearly visible (≥8° angular span)
        if abs(a_end - a_start) < 8.0:
            mid = (a_start + a_end) / 2
            a_start, a_end = mid - 4.0, mid + 4.0
        n = 40
        thetas = list(_linspace(a_start, a_end, n))
        # Opacity floor of 0.55 so weak loci are still visible
        opacity = 0.55 + 0.40 * (count / max_count)
        fillcolor = IG_FILL.get(gene, "#1d4ed8")
        fig.add_trace(go.Scatterpolar(
            theta=thetas + thetas[::-1],
            r=[R_HL_OUT] * n + [R_HL_IN] * n,
            mode="lines",
            fill="toself",
            line=dict(width=0),
            fillcolor=fillcolor,
            opacity=opacity,
            hovertext=f"<b>{gene}</b> — {count} physiological rearrangements (V(D)J / class-switch)",
            hoverinfo="text",
            showlegend=False,
        ))
        # Locus label outside the band
        mid = (a_start + a_end) / 2
        fig.add_trace(go.Scatterpolar(
            theta=[mid], r=[R_HL_OUT + 0.08],
            mode="text",
            text=[f"<b>{gene}</b> <span style='font-size:10px;color:#666'>n={count}</span>"],
            textfont=dict(size=13, color=fillcolor),
            showlegend=False, hoverinfo="skip",
        ))

    # Inter-IG arcs (rare, faint background)
    for c in inter_pairs[:200]:
        a1 = _pos_to_angle(c.chrom_a, c.pos_a, sectors)
        a2 = _pos_to_angle(c.chrom_b, c.pos_b, sectors)
        if a1 is None or a2 is None:
            continue
        thetas, rs = _arc_points(a1, a2, R_IN)
        fig.add_trace(go.Scatterpolar(
            theta=thetas, r=rs,
            mode="lines",
            line=dict(color="#9ecae1", width=1.0),
            opacity=0.35,
            hovertext=(f"IG-IG: {c.gene_a or c.chrom_a}:{c.pos_a} ⇄ "
                       f"{c.gene_b or c.chrom_b}:{c.pos_b}"),
            hoverinfo="text",
            showlegend=False,
        ))

    # ── Canonical translocations — bold red arcs + gene-pair labels ──
    canonical = [c for c in calls if classify_translocation(c) is not None]
    # one arc per unique gene pair per sample (keep best evidence)
    label_set: dict[tuple[str, str], FusionCall] = {}
    for c in canonical:
        key = tuple(sorted([c.gene_a, c.gene_b]))
        cur = label_set.get(key)
        score = (tier_rank.get(c.tier, 3), -(c.split_reads + c.discordant_pairs))
        if cur is None or score < (tier_rank.get(cur.tier, 3),
                                   -(cur.split_reads + cur.discordant_pairs)):
            label_set[key] = c

    # Collect canonical labels by angular bucket so we can stack overlapping ones
    label_bucket: dict[int, list[tuple[float, str]]] = {}
    for (ga, gb), c in label_set.items():
        a1 = _pos_to_angle(c.chrom_a, c.pos_a, sectors)
        a2 = _pos_to_angle(c.chrom_b, c.pos_b, sectors)
        if a1 is None or a2 is None:
            continue
        thetas, rs = _arc_points(a1, a2, R_IN)
        # Heavy red arc
        fig.add_trace(go.Scatterpolar(
            theta=thetas, r=rs,
            mode="lines",
            line=dict(color="#b91c1c", width=3.5),
            opacity=0.9,
            hovertext=(f"<b>{c.gene_a} ↔ {c.gene_b}</b><br>{c.known_partner_source}<br>"
                       f"tier={c.tier} | SR={c.split_reads} PE={c.discordant_pairs}"),
            hoverinfo="text",
            showlegend=False,
        ))
        # Bucket each endpoint label by 6° angular position for overlap stacking
        for theta_deg, gene in ((thetas[0], c.gene_a), (thetas[-1], c.gene_b)):
            bkt = int(theta_deg // 6)
            label_bucket.setdefault(bkt, []).append((theta_deg, gene))

    # Render canonical labels stacked radially when buckets contain >1
    for bkt, items in label_bucket.items():
        # de-duplicate identical genes at the same bucket
        seen: dict[str, float] = {}
        for theta_deg, gene in items:
            if gene not in seen:
                seen[gene] = theta_deg
        for i, (gene, theta_deg) in enumerate(seen.items()):
            r_label = 1.32 + i * 0.06     # stack outwards
            fig.add_trace(go.Scatterpolar(
                theta=[theta_deg], r=[r_label],
                mode="text",
                text=[f"<b>{gene}</b>"],
                textfont=dict(size=14, color="#b91c1c"),
                showlegend=False, hoverinfo="skip",
            ))

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=16, color="#16213e")),
        polar=dict(
            radialaxis=dict(visible=False, range=[0, 1.55]),
            angularaxis=dict(visible=False, direction="clockwise", rotation=90),
            bgcolor="white",
        ),
        height=max(height, 820),
        margin=dict(l=30, r=30, t=60, b=50),
        showlegend=False,
        paper_bgcolor="white",
        annotations=[dict(
            text=("<span style='color:#b91c1c;font-weight:600'>● canonical translocation</span> &nbsp;&nbsp; "
                  "<span style='color:#1d4ed8;font-weight:600'>■ IG locus (V(D)J / class-switch density)</span> &nbsp;&nbsp; "
                  "<span style='color:#7c3aed;font-weight:600'>■ TR locus</span>"),
            x=0.5, y=-0.05, xref="paper", yref="paper",
            showarrow=False, font=dict(size=12),
        )],
    )
    return fig


def _linspace(a: float, b: float, n: int):
    if n <= 1:
        yield a
        return
    step = (b - a) / (n - 1)
    for i in range(n):
        yield a + step * i
