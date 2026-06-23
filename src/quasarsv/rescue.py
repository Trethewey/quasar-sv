"""IG-switch-region rescue.

When a B-cell translocation breaks an IG switch (S-μ, S-γ, S-α), the
chimeric reads can mismap to reference artefact hotspots (e.g. chr2:32916
polyG attractor in GRCh38) instead of their true non-IG partner. The
signature: driver_locus and IG_locus both show SR support to the same
spurious partner. Synthesise the inferred driver-IG fusion from that pair.

On PMBL the polyG attractor absorbs nearly all IGH switch reads. IGL and
IGK can carry more per-gene SR than the truth (IGH) because they are
smaller / more compact, so read-data alone cannot disambiguate. Mitigations:

  * lineage prior — B-cell samples consider only {IGH, IGK, IGL};
    T-cell samples {TRA, TRB, TRG, TRD}.
  * canonical alternatives — when multiple IGs are canonical partners of
    the same driver, emit all of them; lower-scoring ones carry the
    ``ig_partner_ambiguous`` flag.
  * fan-out control — non-top canonical drivers must clear
    ``noncanonical_fanout_ratio × top_canonical_score`` before emitting.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Mapping

from .annotate import (
    GeneEntry, KnownPartner,
    _is_known_pair, load_builtin_loci, load_builtin_partners,
)
from .model import FusionCall


# Lineage-specific IG / TR locus sets
BCELL_IGS = {"IGH", "IGK", "IGL", "IGH_Emu", "IGH_3RR"}
TCELL_IGS = {"TRA", "TRB", "TRG", "TRD"}
IG_LOCI = BCELL_IGS | TCELL_IGS


def _allowed_igs_for_lineage(lineage: str) -> set[str]:
    if lineage == "B":
        return BCELL_IGS
    if lineage == "T":
        return TCELL_IGS
    return BCELL_IGS | TCELL_IGS


@dataclass
class RescueConfig:
    min_artefact_sr_per_side: int = 30
    min_artefact_pe_per_side: int = 0
    tier: str = "T2"
    # Only emit drivers / IGs whose support is ≥ ratio_keep × top in that class
    ratio_keep: float = 0.20
    # Hard cap on emitted pairs per sample (defends against fan-out)
    max_pairs_per_sample: int = 12
    # Default lymphoma lineage assumption (B-cell). Overrideable per-sample.
    lineage: str = "B"
    # Emit additional canonical IG candidates for a driver when artefact-mediated
    # signal cannot disambiguate IGH/IGK/IGL. Marked with ``ig_partner_ambiguous``.
    emit_canonical_alternatives: bool = True
    # Cap on canonical IG alternatives per driver (primary + alternatives)
    max_canonical_igs_per_driver: int = 3
    # Drop low-scoring canonical primaries from non-top drivers to control
    # the artefact-driven fan-out of pseudo-canonical pairs (e.g. BCL2-IGH,
    # CCND1-IGH, MALT1-IGH all firing on the same shared IGH signal when only
    # one driver is actually rearranged). 0.0 disables. Threshold is applied
    # as ``score >= ratio * top_canonical_score``.
    noncanonical_fanout_ratio: float = 0.20
    # Cap on non-canonical pairs (one per driver, capped at this total)
    max_noncanonical_pairs_per_sample: int = 4


def _gene_for_position(loci: list[GeneEntry], chrom: str, pos: int) -> str:
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    for g in loci:
        if g.chrom == chrom and g.start - 50_000 <= pos <= g.end + 50_000:
            return g.gene
    return ""


def _ig_lineage(gene: str) -> str:
    if gene in BCELL_IGS:
        return "B"
    if gene in TCELL_IGS:
        return "T"
    return ""


def rescue_ig_driver_pairs(
    calls: list[FusionCall],
    cfg: RescueConfig | None = None,
    sample_lineage: Mapping[str, str] | None = None,
) -> list[FusionCall]:
    """Append synthetic FusionCall entries to ``calls`` for inferred IG-driver pairs.

    Parameters
    ----------
    calls
        FusionCall list (mutated in place — synthetic entries appended).
    cfg
        :class:`RescueConfig` controlling thresholds, lineage default, and
        canonical-alternative emission.
    sample_lineage
        Optional mapping ``sample_id -> "B" | "T" | "any"``. Overrides
        ``cfg.lineage`` per sample. When absent, ``cfg.lineage`` is used.

    Returns
    -------
    list[FusionCall]
        The same input list, with synthetic rescued calls appended.
    """
    cfg = cfg or RescueConfig()
    sample_lineage = sample_lineage or {}
    loci = load_builtin_loci()
    partners = load_builtin_partners()

    # Group artefact-flagged calls by (sample, non-artefact-side-gene)
    by_sample_gene: dict[tuple[str, str], dict] = defaultdict(
        lambda: {"sr": 0, "pe": 0, "examples": [], "is_ig": False,
                 "is_driver": False, "chrom": "", "pos": 0}
    )

    from .qc import _load_artefact_loci
    art = _load_artefact_loci()

    def in_art(ch: str, p: int) -> bool:
        ch_n = ch[3:] if ch.startswith("chr") else ch
        for chrom, start, end, _ in art:
            if ch_n == chrom and start <= p <= end:
                return True
        return False

    for c in calls:
        if "builtin_artefact_locus" not in c.qc_flags:
            continue
        # Determine which side is the artefact and which is the "real" driver/IG
        a_is_art = in_art(c.chrom_a, c.pos_a)
        b_is_art = in_art(c.chrom_b, c.pos_b)
        if a_is_art and not b_is_art:
            real_chrom, real_pos = c.chrom_b, c.pos_b
            real_gene = c.gene_b
        elif b_is_art and not a_is_art:
            real_chrom, real_pos = c.chrom_a, c.pos_a
            real_gene = c.gene_a
        else:
            continue   # both artefact or neither — nothing to rescue
        if not real_gene:
            real_gene = _gene_for_position(loci, real_chrom, real_pos)
        if not real_gene:
            continue
        key = (c.sample, real_gene)
        d = by_sample_gene[key]
        d["sr"] += c.split_reads
        d["pe"] += c.discordant_pairs
        d["examples"].append(c.fusion_id)
        d["is_ig"] = real_gene in IG_LOCI
        d["is_driver"] = (not d["is_ig"]) and any(
            g.gene == real_gene and g.role == "driver" for g in loci
        )
        d["chrom"], d["pos"] = real_chrom, real_pos

    # Bucket by sample → list of (gene, info)
    by_sample: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for (sample, gene), info in by_sample_gene.items():
        by_sample[sample].append((gene, info))

    new_calls: list[FusionCall] = []
    for sample, items in by_sample.items():
        lineage = sample_lineage.get(sample) or cfg.lineage
        allowed_igs = _allowed_igs_for_lineage(lineage)

        igs_all = [(g, i) for g, i in items
                   if i["is_ig"]
                   and g in allowed_igs
                   and i["sr"] >= cfg.min_artefact_sr_per_side]
        drs_all = [(g, i) for g, i in items
                   if i["is_driver"]
                   and i["sr"] >= cfg.min_artefact_sr_per_side]

        if not igs_all or not drs_all:
            continue

        # Drivers: ratio_keep against the top driver SR
        top_dr = max(i["sr"] for _, i in drs_all)
        drs = [(g, i) for g, i in drs_all
               if i["sr"] >= cfg.ratio_keep * top_dr]

        # IGs: ratio_keep against the top IG SR, BUT keep any IG that is a
        # canonical partner of a kept driver (so the truth survives even when
        # its per-gene SR is dwarfed by an artefact-routed competitor IG).
        top_ig = max(i["sr"] for _, i in igs_all)
        igs_ratio = [(g, i) for g, i in igs_all
                     if i["sr"] >= cfg.ratio_keep * top_ig]
        canonical_keep_names: set[str] = set()
        if cfg.emit_canonical_alternatives:
            for d_gene, _ in drs:
                for ig_gene, _ in igs_all:
                    if _is_known_pair(partners, d_gene, ig_gene):
                        canonical_keep_names.add(ig_gene)
        keep_names = {g for g, _ in igs_ratio} | canonical_keep_names
        igs = [(g, i) for g, i in igs_all if g in keep_names]

        if not igs:
            continue

        # Score candidate pairs; canonical pairs come first, then by min(SR).
        scored: list[tuple[int, float, str, dict, str, dict]] = []
        for d_gene, d_info in drs:
            for ig_gene, ig_info in igs:
                is_kp = _is_known_pair(partners, d_gene, ig_gene) is not None
                # Pair score = min of driver/IG per-gene SR contribution to artefact
                score = float(min(d_info["sr"], ig_info["sr"]))
                scored.append((1 if is_kp else 0, score, d_gene, d_info, ig_gene, ig_info))
        # Sort: canonical first (descending), then score descending
        scored.sort(key=lambda r: (-r[0], -r[1]))

        # Identify the top canonical driver — only this driver emits multiple
        # canonical IG alternatives; other canonical drivers only emit their
        # single top canonical pair, gated on noncanonical_fanout_ratio.
        top_canonical_driver: str | None = None
        top_canonical_score: float = 0.0
        for is_kp, sc, d_gene, *_ in scored:
            if is_kp:
                top_canonical_driver = d_gene
                top_canonical_score = sc
                break

        emitted_for_driver: dict[str, int] = defaultdict(int)
        emitted_pairs: set[tuple[str, str]] = set()
        noncanonical_emitted = 0
        for is_kp_int, score, d_gene, d_info, ig_gene, ig_info in scored:
            if len(new_calls) >= cfg.max_pairs_per_sample:
                break
            if (d_gene, ig_gene) in emitted_pairs:
                continue
            primary_count = emitted_for_driver[d_gene]

            if is_kp_int:
                # Canonical pair
                if primary_count == 0:
                    # First canonical emission for this driver. Apply
                    # fan-out threshold for non-top-canonical drivers so a
                    # single weak shared IG signal does not light up every
                    # driver that happens to be a canonical IG partner.
                    if (top_canonical_driver is not None
                            and d_gene != top_canonical_driver
                            and cfg.noncanonical_fanout_ratio > 0
                            and score < cfg.noncanonical_fanout_ratio * top_canonical_score):
                        continue
                else:
                    # Canonical alt — only the top canonical driver may emit alts.
                    if d_gene != top_canonical_driver:
                        continue
                    if primary_count >= cfg.max_canonical_igs_per_driver:
                        continue
            else:
                # Non-canonical pair: at most one per driver, capped overall.
                if primary_count >= 1:
                    continue
                if noncanonical_emitted >= cfg.max_noncanonical_pairs_per_sample:
                    continue

            fusion_id = f"{sample}__rescue__{d_gene}_{ig_gene}"
            synthetic = FusionCall(
                sample=sample,
                fusion_id=fusion_id,
                chrom_a=_strip(d_info["chrom"]),
                pos_a=d_info["pos"], strand_a=".",
                chrom_b=_strip(ig_info["chrom"]),
                pos_b=ig_info["pos"], strand_b=".",
                sv_type="BND",
                callers_supporting=["forge_scan_rescue"],
                n_callers=1,
                split_reads=min(d_info["sr"], ig_info["sr"]),
                discordant_pairs=min(d_info["pe"], ig_info["pe"]),
                assembly_contigs=0,
                soft_clips=0,
                n_evidence_types=1 + (1 if min(d_info["pe"], ig_info["pe"]) > 0 else 0),
                vaf=0.0,
                precise=False,
                any_pass=True,
                raw_qual_max=float(min(d_info["sr"], ig_info["sr"])),
                gene_a=d_gene, region_a="rearr_hotspot",
                gene_b=ig_gene, region_b="IG_locus",
                in_frame=None,
                known_partner=False,
                known_partner_source="",
                driver_locus=f"{d_gene}-{ig_gene}",
                tier=cfg.tier,
                qc_flags=["inferred_via_artefact_rescue"],
                member_record_ids=d_info["examples"] + ig_info["examples"],
            )
            kp = _is_known_pair(partners, d_gene, ig_gene)
            if kp:
                synthetic.known_partner = True
                synthetic.known_partner_source = f"{kp.source}:{kp.disease}"
                if synthetic.tier == "T2":
                    synthetic.tier = "T1"   # canonical partner + reciprocal artefact signal = T1

            # Ambiguity tag: this driver already had a higher-scoring IG candidate emitted.
            if primary_count >= 1:
                synthetic.qc_flags.append("ig_partner_ambiguous")

            emitted_pairs.add((d_gene, ig_gene))
            emitted_for_driver[d_gene] += 1
            if not is_kp_int:
                noncanonical_emitted += 1
            new_calls.append(synthetic)

    calls.extend(new_calls)
    return calls


def _strip(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c
