"""Pysam-based read-level fusion scanner — split-read + discordant-pair evidence
across a defined set of target loci (default = built-in lymphoma driver loci).

Emits BreakpointCall records compatible with the merger. The "caller" name is
``forge_scan`` so it slots alongside Manta / GRIDSS / Delly / SvABA / FACTERA
in the ensemble.

Design
------
For every read overlapping each target locus:
1. **Split reads.** Read carries an SA (supplementary alignment) tag → take
   the primary side as one breakpoint and the SA alignment as the other.
   Strands derive from the soft-clip side (leading clip = '-' end, trailing
   clip = '+' end).
2. **Discordant pairs.** Read is properly mapped but its mate maps to a
   different chromosome or > `discordant_min_distance` bp away on the same
   chromosome. Pair contributes to the discordant-pair count of a breakpoint
   cluster.

Mate / SA positions are clustered with ±`pos_tolerance` to produce candidate
breakpoint pairs. A cluster supported by ≥`min_split_reads` split reads OR
≥`min_discordant_pairs` discordant pairs is emitted as a BreakpointCall.

Performance
-----------
Targeted-only scanning of ~36 lymphoma loci (~5–10 Mb total) on a 50–60 GB
WGS CRAM completes in minutes via CRAI random access. No whole-genome read
streaming.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

try:
    import pysam
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pysam is required for quasarsv.scanners — install with `pip install quasarsv[bam]` "
        "(may need WSL / Linux for pysam wheels)."
    ) from e

from ..annotate import GeneEntry, load_builtin_loci
from ..model import BreakpointCall, Evidence
from ..parsers.base import normalise_order


@dataclass
class ScannerConfig:
    min_mapq: int = 10                     # hard floor (lowered from 20 — see mapq weighting)
    full_weight_mapq: int = 30             # ≥ this contributes weight 1.0; below tapers linearly
    use_mapq_weighting: bool = True        # GRIDSS-style soft weighting (per docs/precision_techniques.md #10)
    discordant_min_distance: int = 10_000  # initial value; overridden by library_stats when present
    use_adaptive_insert: bool = True       # Delly-style median + 5×MAD threshold (per docs #6)
    pos_tolerance: int = 500
    min_split_reads: int = 2
    min_discordant_pairs: int = 4
    pad_locus_bp: int = 5_000
    max_reads_per_locus: int = 5_000_000   # safety cap
    chr_prefix: bool = True                # CRAM uses chrN names
    # Discard split reads whose clipped segment is adapter read-through or a
    # 2-colour poly-G tail. These align to reference homopolymer attractors and
    # fabricate junction support uniformly across the genome. Off only for
    # diagnosing the artefact channel itself.
    filter_noise_clips: bool = True


@dataclass
class _Cluster:
    chrom_b: str
    pos_b: int
    strand_a: str = "+"
    strand_b: str = "+"
    split_reads: float = 0.0          # float to support MAPQ-weighted contributions
    discordant_pairs: float = 0.0
    pos_a_examples: list[int] = field(default_factory=list)
    mapqs: list[int] = field(default_factory=list)


# Illumina TruSeq adapter and its reverse complement. Read-through into the
# adapter produces a soft-clip that aligns nowhere real.
_ADAPTER = "AGATCGGAAGAGC"
_ADAPTER_RC = "GCTCTTCCGATCT"
# 2-colour chemistry (NovaSeq/NextSeq) encodes "no signal" as G, so low-quality
# read tails degrade into poly-G (poly-C on the reverse strand).
_HOMOPOLYMER_FRAC = 0.70


def _parse_sa_tag(sa: str) -> list[tuple[str, int, str, str]]:
    """Parse EVERY SA entry: 'chrom,pos,strand,CIGAR,mapQ,NM;...'.

    Chimeric reads at IG loci routinely carry more than two alignments, so
    taking only the first entry silently discards real breakpoints.
    """
    out: list[tuple[str, int, str, str]] = []
    for entry in sa.split(";"):
        if not entry:
            continue
        parts = entry.split(",")
        if len(parts) < 4:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        out.append((parts[0], pos, parts[2], parts[3]))
    return out


def _longest_clip_seq(read) -> str:
    """Sequence of the read's longest soft-clipped segment ('' if none)."""
    if not read.cigartuples or read.query_sequence is None:
        return ""
    seq = read.query_sequence
    best = ""
    op, ln = read.cigartuples[0]
    if op == 4 and ln > len(best):
        best = seq[:ln]
    op, ln = read.cigartuples[-1]
    if op == 4 and ln > len(best):
        best = seq[-ln:]
    return best


def is_noise_clip(seq: str) -> bool:
    """True when a soft-clip is a sequencing artefact rather than genomic sequence.

    Adapter read-through and 2-colour poly-G tails align opportunistically to
    reference homopolymer attractors (notably the GRCh38 poly-G at chr2:32,916),
    manufacturing split-read "evidence" that links every locus to that sink at a
    uniform rate — measured at ~220 SR per 10k reads for rearranged drivers and
    unrearranged controls alike. Such clips carry no junction information.
    """
    if not seq:
        return False
    if _ADAPTER[:11] in seq or _ADAPTER_RC[:11] in seq:
        return True
    n = len(seq)
    return (seq.count("G") / n >= _HOMOPOLYMER_FRAC
            or seq.count("C") / n >= _HOMOPOLYMER_FRAC)


def _clip_side(cigar: str) -> str:
    """Return 'L' if leading soft-clip dominates, 'R' if trailing.

    A read with leading clip means the breakpoint joins to upstream sequence
    of the read's mapping; convention: strand '-' on this side. Trailing clip
    → '+' (read continues downstream of breakpoint).
    """
    lead = 0
    trail = 0
    digits = ""
    last_n = 0
    last_op = ""
    for ch in cigar:
        if ch.isdigit():
            digits += ch
        else:
            n = int(digits) if digits else 0
            digits = ""
            if not lead and ch in ("S", "H") and last_n == 0:
                lead = n
            last_n = n
            last_op = ch
    if last_op in ("S", "H"):
        trail = last_n
    if lead >= trail and lead > 0:
        return "L"
    if trail > 0:
        return "R"
    return "."


def _bucket(pos: int, tol: int) -> int:
    return pos // tol


def scan_cram(
    bam_or_cram: str,
    reference_fasta: str,
    sample: str,
    loci: Iterable[GeneEntry] | None = None,
    cfg: ScannerConfig | None = None,
    library_stats_path: str | None = None,
) -> list[BreakpointCall]:
    """Scan a BAM/CRAM and return BreakpointCall list (caller='forge_scan').

    If ``library_stats_path`` is provided, cached library stats are loaded
    from there (or computed + saved if absent) and used to set the adaptive
    ``discordant_min_distance`` per Delly's median + 5 × MAD heuristic
    (see ``library_stats.py``).
    """
    cfg = cfg or ScannerConfig()
    loci = list(loci) if loci is not None else load_builtin_loci()

    # Adaptive insert threshold from library stats (Delly-style)
    if cfg.use_adaptive_insert and library_stats_path:
        from .library_stats import LibraryStats, compute_library_stats
        lib_path = library_stats_path
        if Path(lib_path).exists():
            lib = LibraryStats.load(lib_path)
        else:
            lib = compute_library_stats(bam_or_cram, reference_fasta)
            try:
                lib.save(lib_path)
            except OSError:
                pass    # cache write best-effort
        cfg.discordant_min_distance = lib.discordant_min_distance

    open_kwargs = {}
    if bam_or_cram.endswith(".cram"):
        open_kwargs["reference_filename"] = reference_fasta
        mode = "rc"
    else:
        mode = "rb"
    sam = pysam.AlignmentFile(bam_or_cram, mode, **open_kwargs)

    # CRAM uses 'chr1' style; loci file uses '1' style → adapt
    ref_names = set(sam.references)
    chr_prefix = any(r.startswith("chr") for r in ref_names)

    def to_ref(c: str) -> str:
        if chr_prefix:
            return c if c.startswith("chr") else f"chr{c}"
        return c[3:] if c.startswith("chr") else c

    calls: list[BreakpointCall] = []
    for g in loci:
        ref_chrom = to_ref(g.chrom)
        if ref_chrom not in ref_names:
            continue
        start = max(0, g.start - cfg.pad_locus_bp)
        end = g.end + cfg.pad_locus_bp
        clusters: dict[tuple[str, str, int, str, str], _Cluster] = {}
        n_reads = 0
        for read in sam.fetch(ref_chrom, start, end):
            n_reads += 1
            if n_reads > cfg.max_reads_per_locus:
                break
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < cfg.min_mapq:
                continue

            # MAPQ weighting: contribute fractional support for borderline-mapped
            # reads instead of either dropping (lossy) or counting full (over-
            # counting). Reads at MAPQ ≥ full_weight_mapq contribute 1.0; below
            # tapers linearly to 0 at MAPQ 0.
            if cfg.use_mapq_weighting:
                from .library_stats import mapq_weight
                weight = mapq_weight(read.mapping_quality, cfg.full_weight_mapq)
            else:
                weight = 1.0

            # Split-read evidence via SA tag
            try:
                sa = read.get_tag("SA")
            except KeyError:
                sa = ""
            if sa and cfg.filter_noise_clips and is_noise_clip(_longest_clip_seq(read)):
                # Adapter read-through / poly-G tail: not a junction. Drop the
                # split-read evidence, and drop the read entirely rather than
                # letting it fall through into the discordant path, where the
                # same artefactual mapping would be recounted as a pair.
                sa = ""
                continue
            if sa:
                sa_entries = _parse_sa_tag(sa)
                primary_clip = _clip_side(read.cigarstring or "")
                # An unclipped read cannot tell us which side the junction is
                # on. Abstain ('.') rather than defaulting to '+', which would
                # report a fabricated orientation as a measured one.
                strand_a = {"L": "-", "R": "+"}.get(primary_clip, ".")
                if primary_clip == "L":
                    pa = read.reference_start
                else:
                    pa = read.reference_end or read.reference_start
                for sa_chrom, sa_pos, sa_strand, sa_cigar in sa_entries:
                    sa_clip = _clip_side(sa_cigar)
                    strand_b = {"L": "-", "R": "+"}.get(sa_clip, ".")
                    pb = sa_pos
                    key = (sa_chrom, "SR", _bucket(pb, cfg.pos_tolerance), strand_a, strand_b)
                    cl = clusters.get(key)
                    if cl is None:
                        cl = _Cluster(chrom_b=sa_chrom, pos_b=pb,
                                      strand_a=strand_a, strand_b=strand_b)
                        clusters[key] = cl
                    cl.split_reads += weight
                    cl.pos_a_examples.append(pa)
                    cl.mapqs.append(read.mapping_quality)
                if sa_entries:
                    continue   # don't double-count this read as discordant too

            # Discordant-pair evidence
            if read.is_paired and not read.mate_is_unmapped:
                mate_chrom = sam.get_reference_name(read.next_reference_id) or ""
                same_chrom = (read.reference_id == read.next_reference_id)
                far_enough = (
                    not same_chrom
                    or abs(read.next_reference_start - read.reference_start) >= cfg.discordant_min_distance
                )
                if far_enough:
                    strand_a = "-" if read.is_reverse else "+"
                    strand_b = "-" if read.mate_is_reverse else "+"
                    pa = read.reference_end or read.reference_start
                    pb = read.next_reference_start
                    key = (mate_chrom, "PE", _bucket(pb, cfg.pos_tolerance), strand_a, strand_b)
                    cl = clusters.get(key)
                    if cl is None:
                        cl = _Cluster(chrom_b=mate_chrom, pos_b=pb,
                                      strand_a=strand_a, strand_b=strand_b)
                        clusters[key] = cl
                    cl.discordant_pairs += weight
                    cl.pos_a_examples.append(pa)
                    cl.mapqs.append(read.mapping_quality)

        for key, cl in clusters.items():
            if (cl.split_reads < cfg.min_split_reads
                    and cl.discordant_pairs < cfg.min_discordant_pairs):
                continue
            # representative chrom_a = locus chromosome; pos_a = median of supporting reads
            cl.pos_a_examples.sort()
            pa_rep = cl.pos_a_examples[len(cl.pos_a_examples) // 2]
            chrom_a = ref_chrom[3:] if ref_chrom.startswith("chr") else ref_chrom
            chrom_b = cl.chrom_b[3:] if cl.chrom_b.startswith("chr") else cl.chrom_b
            ca, pa, sa_strand, cb, pb, sb_strand = normalise_order(
                chrom_a, pa_rep, cl.strand_a, chrom_b, cl.pos_b, cl.strand_b
            )
            # Round MAPQ-weighted floats to ints for the downstream schema.
            sr_int = max(0, int(round(cl.split_reads)))
            pe_int = max(0, int(round(cl.discordant_pairs)))
            ev = Evidence(
                caller="forge_scan",
                split_reads=sr_int,
                discordant_pairs=pe_int,
                assembly_contigs=0,
                soft_clips=sr_int,
                mapq=int(sum(cl.mapqs) / len(cl.mapqs)) if cl.mapqs else 0,
                vaf=0.0,
                filter_pass=True,
                precise=cl.split_reads > 0,
                raw_qual=float(max(cl.split_reads, cl.discordant_pairs)),
            )
            calls.append(BreakpointCall(
                sample=sample,
                chrom_a=ca, pos_a=pa, strand_a=sa_strand,
                chrom_b=cb, pos_b=pb, strand_b=sb_strand,
                sv_type="BND",
                evidence=ev,
                record_id=f"forge_scan_{g.gene}_{ca}_{pa}_{cb}_{pb}",
            ))
    sam.close()
    return calls


def scan_to_breakpoint_calls(
    bam_or_cram: str,
    reference_fasta: str,
    sample: str,
    output_tsv: str | None = None,
    loci: Iterable[GeneEntry] | None = None,
    cfg: ScannerConfig | None = None,
) -> list[BreakpointCall]:
    """Convenience: scan + optionally write a per-caller TSV summary."""
    calls = scan_cram(bam_or_cram, reference_fasta, sample, loci=loci, cfg=cfg)
    if output_tsv:
        import csv
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sample", "caller", "chrom_a", "pos_a", "strand_a",
                        "chrom_b", "pos_b", "strand_b",
                        "split_reads", "discordant_pairs", "mapq", "record_id"])
            for b in calls:
                e = b.evidence
                w.writerow([b.sample, e.caller, b.chrom_a, b.pos_a, b.strand_a,
                            b.chrom_b, b.pos_b, b.strand_b,
                            e.split_reads, e.discordant_pairs, e.mapq, b.record_id])
    return calls
