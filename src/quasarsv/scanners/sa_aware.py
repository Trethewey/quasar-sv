"""SA-tag aware scanner for reads landing in a reference artefact locus.

The premise here needs stating carefully, because it was originally wrong. The
GRCh38 poly-G attractor at chr2:32,916 does NOT preferentially absorb chimeric
IG-switch fragments. It absorbs **2-colour-chemistry poly-G tails and adapter
read-through from the entire library**, at a uniform rate: every locus measured
on the WGS validation cohort sheds ~200-280 such reads per 10k, rearranged or
not (``quasar_development/accuracy_audit_2026-07-16/artefact_specificity.py``).

So harvesting SA tags from reads inside the attractor does not "reveal the real
translocation partners" — for poly-G reads it merely enumerates whichever locus
each junk read originated from, i.e. the whole panel. That is the source of the
artefact's apparent promiscuity.

Note the read geometry: at the artefact the ALIGNED segment is the poly-G run
and the soft-clip holds the real genomic sequence — the mirror image of the
driver-side view. So the clip-based noise filter used by ``cram_scanner`` does
not apply here; the aligned segment must be tested instead.

What survives the filter is a genuine chimera that happens to overlap the
artefact window. Those are still emitted, as ordinary candidate breakpoints.

Output is a list of BreakpointCall records:
  chrom_a = partner chromosome (from SA tag)
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
from .cram_scanner import _parse_sa_tag, _clip_side, is_noise_clip


@dataclass
class SAScannerConfig:
    min_mapq: int = 20
    pos_tolerance: int = 500
    min_split_reads: int = 5
    pad_locus_bp: int = 200
    # Polymorphic artefact loci attract millions of reads — cap aggressively
    # since even 50k reads typically yields rich SA-tag clustering signal.
    max_reads_per_locus: int = 200_000
    # Reject reads whose ALIGNED segment inside the artefact is a poly-G/C run
    # or adapter. Those are sequencing junk from an arbitrary locus, and their
    # SA tag names that origin locus rather than any translocation partner.
    filter_noise_alignments: bool = True


def _aligned_seq(read) -> str:
    """The portion of the read actually aligned here (soft-clips excluded)."""
    if read.query_sequence is None or not read.cigartuples:
        return ""
    start = 0
    op, ln = read.cigartuples[0]
    if op == 4:
        start = ln
    end = len(read.query_sequence)
    op, ln = read.cigartuples[-1]
    if op == 4:
        end -= ln
    return read.query_sequence[start:end] if end > start else ""


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
            # The read is only informative if what aligned HERE is real sequence.
            # A poly-G run or adapter aligned into the attractor tells us nothing
            # except which locus the junk read came from.
            if cfg.filter_noise_alignments and is_noise_clip(_aligned_seq(read)):
                continue
            primary_clip = _clip_side(read.cigarstring or "")
            strand_b = "-" if primary_clip == "L" else "+"
            for sa_chrom, sa_pos, sa_strand, sa_cigar in _parse_sa_tag(sa):
                # Don't re-cluster intra-artefact
                if sa_chrom == ref_chrom and s - 1000 <= sa_pos <= e + 1000:
                    continue
                sa_clip = _clip_side(sa_cigar)
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
