"""Canonical translocation classes for grouping reports.

Each class defines a cytogenetic / disease-named bucket (e.g. t(14;18) IGH-BCL2)
plus the gene pairs that fall into it. Driving the brochure and cohort report's
"one table per translocation" layout.
"""
from __future__ import annotations

from dataclasses import dataclass

from .model import FusionCall


@dataclass
class TranslocationClass:
    key: str                       # short id, e.g. 't_14_18_IGH_BCL2'
    cytoband: str                  # 't(14;18)(q32;q21)'
    label: str                     # 'IGH-BCL2 — Follicular / DH-DLBCL'
    disease: str                   # 'FL, DLBCL (DH/TH)'
    gene_pairs: list[tuple[str, str]]   # canonical gene-pair list


CANONICAL_CLASSES: list[TranslocationClass] = [
    TranslocationClass(
        "t_8_14_MYC", "t(8;14)(q24;q32) and variants",
        "MYC ↔ IGH/IGK/IGL — Burkitt / DH-DLBCL",
        "Burkitt lymphoma; DH/TH-DLBCL",
        [("IGH", "MYC"), ("IGK", "MYC"), ("IGL", "MYC")],
    ),
    TranslocationClass(
        "t_14_18_BCL2", "t(14;18)(q32;q21)",
        "IGH ↔ BCL2 — Follicular / DH-DLBCL",
        "Follicular lymphoma; DH/TH-DLBCL",
        [("IGH", "BCL2")],
    ),
    TranslocationClass(
        "t_3_14_BCL6", "t(3;14)(q27;q32) and variants",
        "BCL6 ↔ IGH/IGK/IGL — DLBCL (incl. PMBL)",
        "Diffuse large B-cell lymphoma",
        [("IGH", "BCL6"), ("IGK", "BCL6"), ("IGL", "BCL6")],
    ),
    TranslocationClass(
        "t_11_14_CCND1", "t(11;14)(q13;q32)",
        "IGH ↔ CCND1 — Mantle cell lymphoma",
        "Mantle cell lymphoma",
        [("IGH", "CCND1")],
    ),
    TranslocationClass(
        "t_6_14_CCND3", "t(6;14)(p21;q32)",
        "IGH ↔ CCND3 — MCL variant",
        "MCL variant",
        [("IGH", "CCND3")],
    ),
    TranslocationClass(
        "MALT", "t(11;18) and variants",
        "MALT1 / BCL10 / API2 rearrangements",
        "MALT lymphoma",
        [("IGH", "MALT1"), ("API2", "MALT1"), ("IGH", "BCL10")],
    ),
    TranslocationClass(
        "ALK_fusions", "t(2;5) and variants",
        "ALK fusions (NPM1-ALK, TPM3-ALK, …)",
        "Anaplastic large-cell lymphoma",
        [("NPM1", "ALK"), ("TPM3", "ALK"), ("TPM4", "ALK"),
         ("ATIC", "ALK"), ("TFG", "ALK")],
    ),
    TranslocationClass(
        "DUSP22_IRF4", "t(6;6)(p25;p25)",
        "DUSP22 ↔ IRF4 — ALK-neg ALCL / cutaneous TCL",
        "ALK-negative ALCL; cutaneous T-cell lymphoma",
        [("DUSP22", "IRF4")],
    ),
    TranslocationClass(
        "t_6_14_IRF4", "t(6;14)(p25;q32)",
        "IGH ↔ IRF4 — myeloma / rare DLBCL",
        "Myeloma; rare DLBCL",
        [("IGH", "IRF4")],
    ),
    TranslocationClass(
        "t_14_19_BCL3", "t(14;19)(q32;q13)",
        "IGH ↔ BCL3 — CLL",
        "Chronic lymphocytic leukaemia",
        [("IGH", "BCL3")],
    ),
    TranslocationClass(
        "t_9_14_PAX5", "t(9;14)(p13;q32)",
        "PAX5 ↔ IGH — LPL / some DLBCL",
        "LPL; some DLBCL",
        [("IGH", "PAX5")],
    ),
    TranslocationClass(
        "BCL6_promiscuous", "BCL6 non-IG partners",
        "BCL6 ↔ non-IG promiscuous partners",
        "DLBCL (BCL6 promiscuous)",
        [("BCL6", "TBL1XR1"), ("BCL6", "LCP1"), ("NFKBIE", "BCL6")],
    ),
]


def _pair_key(g1: str, g2: str) -> frozenset:
    return frozenset((g1, g2))


_PAIR_LOOKUP: dict[frozenset, TranslocationClass] = {}
for _cls in CANONICAL_CLASSES:
    for a, b in _cls.gene_pairs:
        _PAIR_LOOKUP[_pair_key(a, b)] = _cls


def classify_translocation(call: FusionCall) -> TranslocationClass | None:
    """Return the canonical translocation class matching `call.gene_a/gene_b`."""
    if not (call.gene_a and call.gene_b):
        return None
    return _PAIR_LOOKUP.get(_pair_key(call.gene_a, call.gene_b))


def group_by_translocation(calls: list[FusionCall]) -> dict[str, list[FusionCall]]:
    """Group canonical calls by translocation class key. Non-canonical calls
    fall under 'OTHER_canonical' if both genes are known driver/IG but not in
    a defined class. Non-canonical, non-IG-driver calls are NOT included here."""
    out: dict[str, list[FusionCall]] = {}
    for c in calls:
        cls = classify_translocation(c)
        if cls is not None:
            out.setdefault(cls.key, []).append(c)
    return out


def is_ig_involved(call: FusionCall) -> bool:
    """Either breakpoint in an IG or TR locus."""
    from .classify import IG_TR_LOCI
    return (call.gene_a in IG_TR_LOCI) or (call.gene_b in IG_TR_LOCI)
