"""Unified data model used across the pipeline.

Caller parsers produce ``BreakpointCall``. The merger groups them into
``FusionCandidate`` clusters. Annotation enriches each candidate into a
``FusionCall`` with gene-level context. Plots and reports consume ``FusionCall``.

No silent coercions. Field names match the TSV schema exactly.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional
import csv
import json


# Evidence types — these are independent, not just per-caller votes.
EVIDENCE_TYPES = ("split_read", "discordant_pair", "assembly_contig", "soft_clip")


@dataclass
class Evidence:
    """Per-caller per-breakpoint raw evidence as extracted from a caller VCF."""

    caller: str                     # 'manta' | 'gridss' | 'delly' | 'svaba' | 'factera'
    split_reads: int = 0            # SR / ASSR / DV
    discordant_pairs: int = 0       # PE / ASRP / DR
    assembly_contigs: int = 0       # ASC / IC / assembled contigs supporting BP
    soft_clips: int = 0             # BSC / soft-clip count if reported
    mapq: int = 0                   # max / median MAPQ reported by caller
    vaf: float = 0.0                # caller-reported VAF if any
    filter_pass: bool = False       # FILTER == PASS
    precise: bool = True            # ! IMPRECISE flag
    raw_qual: float = 0.0           # QUAL field


@dataclass
class BreakpointCall:
    """A single caller's call at one breakpoint pair (BND, INV, DUP, DEL).

    `chrom_a/pos_a/strand_a` is by convention the lower-sorted breakpoint.
    Strand convention: '+' = breakpoint joins to sequence downstream of pos;
                       '-' = breakpoint joins to sequence upstream of pos.
    SVTYPE BND uses ALT-bracket notation to derive strands.
    """

    sample: str
    chrom_a: str
    pos_a: int
    strand_a: str
    chrom_b: str
    pos_b: int
    strand_b: str
    sv_type: str                    # BND | DEL | DUP | INV | INS | TRA
    evidence: Evidence
    cipos_a: tuple[int, int] = (0, 0)  # confidence interval around pos_a
    cipos_b: tuple[int, int] = (0, 0)  # confidence interval around pos_b
    record_id: str = ""             # original VCF ID for traceability


@dataclass
class FusionCandidate:
    """Cluster of `BreakpointCall`s that the merger considers the same event."""

    sample: str
    chrom_a: str
    pos_a: int
    strand_a: str
    chrom_b: str
    pos_b: int
    strand_b: str
    sv_type: str
    evidences: list[Evidence] = field(default_factory=list)
    member_callers: list[str] = field(default_factory=list)
    member_record_ids: list[str] = field(default_factory=list)

    def evidence_summary(self) -> dict:
        """Sum independent evidence types across all member callers (de-duped by type).

        Two callers seeing the same split-read cluster is one piece of evidence.
        A caller's split-reads + another's assembly contig is two pieces.
        We take the MAX per evidence-type as the representative count.
        """
        summary = {t: 0 for t in EVIDENCE_TYPES}
        for ev in self.evidences:
            summary["split_read"] = max(summary["split_read"], ev.split_reads)
            summary["discordant_pair"] = max(summary["discordant_pair"], ev.discordant_pairs)
            summary["assembly_contig"] = max(summary["assembly_contig"], ev.assembly_contigs)
            summary["soft_clip"] = max(summary["soft_clip"], ev.soft_clips)
        return summary

    def n_independent_evidence_types(self) -> int:
        s = self.evidence_summary()
        return sum(1 for v in s.values() if v > 0)

    def max_vaf(self) -> float:
        return max((e.vaf for e in self.evidences), default=0.0)

    def any_precise(self) -> bool:
        return any(e.precise for e in self.evidences)

    def any_filter_pass(self) -> bool:
        return any(e.filter_pass for e in self.evidences)


@dataclass
class FusionCall:
    """Annotated, tiered fusion event — the unit reports and plots consume."""

    sample: str
    fusion_id: str                  # stable string id per sample
    chrom_a: str
    pos_a: int
    strand_a: str
    chrom_b: str
    pos_b: int
    strand_b: str
    sv_type: str

    # Per-caller / evidence
    callers_supporting: list[str] = field(default_factory=list)
    n_callers: int = 0
    split_reads: int = 0
    discordant_pairs: int = 0
    assembly_contigs: int = 0
    soft_clips: int = 0
    n_evidence_types: int = 0
    vaf: float = 0.0
    precise: bool = False
    any_pass: bool = False
    raw_qual_max: float = 0.0

    # Annotation
    gene_a: str = ""
    gene_b: str = ""
    region_a: str = ""              # exon_N | intron_N | UTR | upstream | downstream
    region_b: str = ""
    in_frame: Optional[bool] = None
    known_partner: bool = False     # canonical lymphoma fusion partner pair
    known_partner_source: str = ""  # COSMIC | Mitelman | lymphoma_panel
    driver_locus: str = ""          # IGH | MYC | BCL2 | BCL6 | CCND1 | ...

    # Tiering and QC flags
    tier: str = "T3"                # T1 | T2 | T3 | filtered
    qc_flags: list[str] = field(default_factory=list)
    # Event classification — physiological vs somatic
    event_class: str = ""           # IG_intra | IG_IG | IG_driver_canonical |
                                    # IG_driver_novel | IG_intergenic |
                                    # driver_driver | driver_intra |
                                    # driver_intergenic | intergenic

    # Provenance
    member_record_ids: list[str] = field(default_factory=list)


# ---- TSV I/O ----

FUSION_CALL_COLUMNS = [
    "sample", "fusion_id", "tier", "event_class",
    "chrom_a", "pos_a", "strand_a", "gene_a", "region_a",
    "chrom_b", "pos_b", "strand_b", "gene_b", "region_b",
    "sv_type", "driver_locus", "known_partner", "known_partner_source",
    "in_frame",
    "n_callers", "callers_supporting",
    "n_evidence_types", "split_reads", "discordant_pairs",
    "assembly_contigs", "soft_clips",
    "vaf", "precise", "any_pass", "raw_qual_max",
    "qc_flags", "member_record_ids",
]


def _scalar(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)


def write_fusion_calls_tsv(calls: list[FusionCall], path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(FUSION_CALL_COLUMNS)
        for c in calls:
            d = asdict(c)
            w.writerow([_scalar(d.get(col, "")) for col in FUSION_CALL_COLUMNS])


def read_fusion_calls_tsv(path: str) -> list[FusionCall]:
    out: list[FusionCall] = []
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            kwargs = dict(row)
            for k in ("pos_a", "pos_b", "n_callers", "n_evidence_types",
                      "split_reads", "discordant_pairs", "assembly_contigs", "soft_clips"):
                kwargs[k] = int(float(kwargs.get(k) or 0))
            for k in ("vaf", "raw_qual_max"):
                kwargs[k] = float(kwargs.get(k) or 0.0)
            for k in ("precise", "any_pass", "known_partner"):
                kwargs[k] = (kwargs.get(k, "").lower() == "true")
            f = kwargs.get("in_frame", "")
            kwargs["in_frame"] = (None if f == "" else f.lower() == "true")
            for k in ("callers_supporting", "qc_flags", "member_record_ids"):
                v = kwargs.get(k) or ""
                kwargs[k] = [s for s in v.split(",") if s]
            out.append(FusionCall(**kwargs))  # type: ignore[arg-type]
    return out


def write_fusion_calls_json(calls: list[FusionCall], path: str) -> None:
    """JSON sibling for downstream tools that prefer structured payloads."""
    payload = [asdict(c) for c in calls]
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
