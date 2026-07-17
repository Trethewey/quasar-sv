"""Chromosome-level SA-tag partner inference — REMOVED, it fabricated coordinates.

The idea was: the position-clustered SA scanner misses diffuse partner signals,
because IG V regions span >1 Mb and SA-tag positions spread too thinly to reach
a per-position cluster threshold. So aggregate at chromosome level instead, and
emit an ``artefact_locus <-> partner`` call for any chromosome carrying >=2% of
the artefact's SA tags.

Two independent defects make that unsalvageable:

1. **The emitted breakpoint is invented.** The call's position was the MEDIAN of
   every SA position on that chromosome. Those positions are unrelated to each
   other, so their median is an arbitrary coordinate that need not sit near any
   real junction — and ``_nearest_gene`` then named whichever driver or IG locus
   happened to lie closest to it. A fabricated coordinate that lands near IGH
   yields a confident-looking IGH partner.
2. **The input is not translocation signal.** Reads inside the chr2:32,916
   attractor are 2-colour poly-G tails and adapter read-through drawn uniformly
   from the whole library (measured: every locus sheds ~200-280 per 10k reads,
   rearranged or not). Ranking their SA chromosomes therefore ranks chromosomes
   roughly by size, so every large chromosome clears a 2% share.

A real diffuse-partner signal is recovered from discordant mates, which map
independently of the junction sequence and carry a true coordinate. That is what
``cram_scanner``'s discordant path does.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

try:
    import pysam
except ImportError as e:  # pragma: no cover
    raise ImportError("pysam required") from e

from ..annotate import GeneEntry, load_builtin_loci
from ..model import BreakpointCall, Evidence
from ..parsers.base import normalise_order
from ..qc import _load_artefact_loci
from .cram_scanner import _parse_sa_tag, _clip_side


@dataclass
class ChromInferenceConfig:
    min_mapq: int = 20
    max_reads_per_artefact: int = 200_000
    min_fraction_of_sa: float = 0.02     # ≥2 % of all SA tags
    min_sa_count: int = 50               # absolute floor


# Loci of clinical interest for partner identification
# (chromosome name → list of (gene, start, end))
def _index_driver_chroms(loci: list[GeneEntry]) -> dict[str, list[GeneEntry]]:
    idx: dict[str, list[GeneEntry]] = {}
    for g in loci:
        idx.setdefault(g.chrom, []).append(g)
    return idx


def _nearest_gene(idx, chrom: str, pos: int, pad: int = 200_000):
    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    for g in idx.get(chrom, []):
        if g.start - pad <= pos <= g.end + pad:
            return g.gene
    return ""


def scan_artefacts_chrom_inference(*_args, **_kwargs) -> list[BreakpointCall]:
    """Removed: emitted a median-of-unrelated-positions as a breakpoint.

    See the module docstring. Diffuse partner signal comes from discordant
    mates via ``cram_scanner``, which carry real coordinates.
    """
    raise NotImplementedError(
        "scan_artefacts_chrom_inference has been removed: it emitted the MEDIAN "
        "of unrelated SA positions on a chromosome as a breakpoint coordinate, "
        "then named the nearest gene as the partner. Its input (reads inside the "
        "chr2:32,916 poly-G attractor) is library-wide adapter/poly-G noise, not "
        "translocation signal. Use the discordant-pair path in cram_scanner."
    )


def _removed_scan_artefacts_chrom_inference(
    bam_or_cram: str,
    reference_fasta: str,
    sample: str,
    cfg: ChromInferenceConfig | None = None,
) -> list[BreakpointCall]:
    cfg = cfg or ChromInferenceConfig()
    art = _load_artefact_loci()
    if not art:
        return []
    loci = load_builtin_loci()
    drv_idx = _index_driver_chroms(loci)

    open_kwargs = {}
    mode = "rc" if bam_or_cram.endswith(".cram") else "rb"
    if mode == "rc":
        open_kwargs["reference_filename"] = reference_fasta
    sam = pysam.AlignmentFile(bam_or_cram, mode, **open_kwargs)
    ref_names = set(sam.references)
    chr_prefix = any(r.startswith("chr") for r in ref_names)
    def to_ref(c: str) -> str:
        return (c if c.startswith("chr") else f"chr{c}") if chr_prefix \
               else (c[3:] if c.startswith("chr") else c)

    calls: list[BreakpointCall] = []
    for art_chrom, art_start, art_end, _notes in art:
        ref_chrom = to_ref(art_chrom)
        if ref_chrom not in ref_names:
            continue
        per_chrom_sa: dict[str, list[int]] = {}
        n_reads = 0
        for read in sam.fetch(ref_chrom, art_start, art_end):
            n_reads += 1
            if n_reads > cfg.max_reads_per_artefact:
                break
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < cfg.min_mapq:
                continue
            try:
                sa = read.get_tag("SA")
            except KeyError:
                continue
            parsed = _parse_sa_tag(sa)
            if parsed is None:
                continue
            sa_chrom, sa_pos, _, _ = parsed
            per_chrom_sa.setdefault(sa_chrom, []).append(sa_pos)

        total_sa = sum(len(v) for v in per_chrom_sa.values())
        if total_sa == 0:
            continue
        # Exclude intra-artefact SAs from the partner ranking
        scored = sorted(
            ((c, len(v), v) for c, v in per_chrom_sa.items() if c != ref_chrom),
            key=lambda r: -r[1],
        )
        for chrom_b, n_sa, positions in scored:
            if n_sa < cfg.min_sa_count:
                continue
            if n_sa / total_sa < cfg.min_fraction_of_sa:
                continue
            # Representative position = median
            positions.sort()
            pos_b = positions[len(positions) // 2]
            gene = _nearest_gene(drv_idx, chrom_b, pos_b)
            cb = chrom_b[3:] if chrom_b.startswith("chr") else chrom_b
            ca = ref_chrom[3:] if ref_chrom.startswith("chr") else ref_chrom
            chrom_a, pos_a, sa_strand, chrom_b_n, pos_b_n, sb_strand = normalise_order(
                ca, art_start + (art_end - art_start) // 2, "+",
                cb, pos_b, "+",
            )
            ev = Evidence(
                caller="quasar_chrom_sa",
                split_reads=n_sa,
                discordant_pairs=0,
                assembly_contigs=0,
                soft_clips=n_sa,
                mapq=0,
                vaf=0.0,
                filter_pass=True,
                precise=False,                       # chrom-level only
                raw_qual=float(n_sa),
            )
            calls.append(BreakpointCall(
                sample=sample,
                chrom_a=chrom_a, pos_a=pos_a, strand_a=sa_strand,
                chrom_b=chrom_b_n, pos_b=pos_b_n, strand_b=sb_strand,
                sv_type="BND", evidence=ev,
                record_id=f"quasar_chrom_sa_{chrom_a}_{chrom_b_n}_{gene or 'NA'}",
            ))
    sam.close()
    return calls
