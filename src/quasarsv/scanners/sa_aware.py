"""SA-tag aware scanner — extracts the TRUE partner of artefact-region reads.

Polymorphic mapping artefacts (e.g. chr2:32916xxx polyG attractor on GRCh38)
absorb soft-clipped fragments of chimeric reads whose real partner is an IG
switch region or a driver oncogene. The SA tag of any read mapping into the
artefact records where the OTHER end of that chimera actually mapped. Scanning
artefact loci and harvesting SA-tag positions therefore reveals the real
translocation partners directly, without the heuristic pairing in `rescue.py`.

Output is again a list of BreakpointCall records:
  chrom_a = real partner chromosome (from SA tag)
  chrom_b = artefact chromosome
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

try:
    import pysam
except ImportError as e:  # pragma: no cover
    raise ImportError("pysam required") from e

from ..model import BreakpointCall, Evidence
from ..parsers.base import normalise_order
from ..qc import _load_artefact_loci
from .cram_scanner import _parse_sa_tag, _clip_side


@dataclass
class SAScannerConfig:
    min_mapq: int = 20
    pos_tolerance: int = 500
    min_split_reads: int = 5
    pad_locus_bp: int = 200
    # Polymorphic artefact loci attract millions of reads — cap aggressively
    # since even 50k reads typically yields rich SA-tag clustering signal.
    max_reads_per_locus: int = 200_000


def scan_artefacts_sa(
    bam_or_cram: str,
    reference_fasta: str,
    sample: str,
    cfg: SAScannerConfig | None = None,
) -> list[BreakpointCall]:
    """Scan the built-in artefact loci and cluster the SA-tag landing sites
    of every read found there. Each strong cluster is emitted as a candidate
    `artefact ↔ real_partner` BreakpointCall, caller=`forge_scan_sa`."""
    cfg = cfg or SAScannerConfig()
    art = _load_artefact_loci()
    if not art:
        return []

    open_kwargs = {}
    mode = "rc" if bam_or_cram.endswith(".cram") else "rb"
    if mode == "rc":
        open_kwargs["reference_filename"] = reference_fasta
    sam = pysam.AlignmentFile(bam_or_cram, mode, **open_kwargs)

    ref_names = set(sam.references)
    chr_prefix = any(r.startswith("chr") for r in ref_names)
    def to_ref(c: str) -> str:
        if chr_prefix:
            return c if c.startswith("chr") else f"chr{c}"
        return c[3:] if c.startswith("chr") else c

    out: list[BreakpointCall] = []
    for chrom, start, end, _notes in art:
        ref_chrom = to_ref(chrom)
        if ref_chrom not in ref_names:
            continue
        s = max(0, start - cfg.pad_locus_bp)
        e = end + cfg.pad_locus_bp
        # Cluster by (real_chrom, bucket)
        clusters: dict[tuple[str, int, str, str], dict] = {}
        n_reads = 0
        for read in sam.fetch(ref_chrom, s, e):
            n_reads += 1
            if n_reads > cfg.max_reads_per_locus:
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
            sa_chrom, sa_pos, sa_strand, sa_cigar = parsed
            # Don't re-cluster intra-artefact
            if sa_chrom == ref_chrom and s - 1000 <= sa_pos <= e + 1000:
                continue
            primary_clip = _clip_side(read.cigarstring or "")
            sa_clip = _clip_side(sa_cigar)
            strand_b = "-" if primary_clip == "L" else "+"
            strand_a = "-" if sa_clip == "L" else "+"
            key = (sa_chrom, sa_pos // cfg.pos_tolerance, strand_a, strand_b)
            d = clusters.setdefault(key, dict(
                sa_chrom=sa_chrom, sa_pos=sa_pos,
                strand_a=strand_a, strand_b=strand_b,
                sr=0, mapqs=[], art_positions=[]))
            d["sr"] += 1
            d["mapqs"].append(read.mapping_quality)
            d["art_positions"].append(read.reference_start)

        for key, d in clusters.items():
            if d["sr"] < cfg.min_split_reads:
                continue
            chrom_a = d["sa_chrom"][3:] if d["sa_chrom"].startswith("chr") else d["sa_chrom"]
            chrom_b = ref_chrom[3:] if ref_chrom.startswith("chr") else ref_chrom
            ca, pa, sa_strand_n, cb, pb, sb_strand_n = normalise_order(
                chrom_a, d["sa_pos"], d["strand_a"],
                chrom_b, d["art_positions"][len(d["art_positions"]) // 2], d["strand_b"],
            )
            ev = Evidence(
                caller="forge_scan_sa",
                split_reads=d["sr"],
                discordant_pairs=0,
                assembly_contigs=0,
                soft_clips=d["sr"],
                mapq=int(sum(d["mapqs"]) / len(d["mapqs"])) if d["mapqs"] else 0,
                vaf=0.0,
                filter_pass=True,
                precise=True,
                raw_qual=float(d["sr"]),
            )
            out.append(BreakpointCall(
                sample=sample, chrom_a=ca, pos_a=pa, strand_a=sa_strand_n,
                chrom_b=cb, pos_b=pb, strand_b=sb_strand_n,
                sv_type="BND", evidence=ev,
                record_id=f"forge_scan_sa_{cb}_{pb}_{ca}_{pa}",
            ))
    sam.close()
    return out
