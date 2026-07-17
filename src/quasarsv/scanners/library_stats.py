"""Per-library statistics computed once per BAM/CRAM at the start of a scan.

Cheap pass that samples 100k reads and infers:

* insert-size median + median-absolute-deviation (MAD)
* soft-clip length distribution
* MAPQ distribution

These feed the adaptive ``discordant_min_distance`` (Delly-style "5 × MAD"
heuristic) and the MAPQ-as-weight scoring (GRIDSS-style — see
``cram_scanner.py``).

All stats are JSON-serialisable so we can cache them per sample under
``output/<sample>/library_stats.json`` and avoid the pass on resume.
"""
from __future__ import annotations

import json
import math
import statistics
from dataclasses import asdict, dataclass, field
from pathlib import Path

try:
    import pysam
except ImportError as e:  # pragma: no cover
    raise ImportError("pysam required for library_stats") from e


@dataclass
class LibraryStats:
    """Per-BAM library inference."""

    n_reads_sampled: int = 0
    insert_median: float = 0.0
    insert_mad: float = 0.0
    insert_p95: float = 0.0
    insert_p99: float = 0.0
    softclip_median: float = 0.0
    softclip_p95: float = 0.0
    mapq_median: float = 0.0
    mapq_p10: float = 0.0
    pct_supplementary: float = 0.0
    pct_duplicate: float = 0.0

    @property
    def discordant_min_distance(self) -> int:
        """The Delly-style adaptive threshold: median + 5 × MAD.

        Falls back to 10_000 when no insert sizes were sampled.
        """
        if self.insert_median <= 0:
            return 10_000
        return int(self.insert_median + 5 * self.insert_mad)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

    @classmethod
    def from_json(cls, s: str) -> "LibraryStats":
        return cls(**json.loads(s))

    @classmethod
    def load(cls, path: str | Path) -> "LibraryStats":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    def save(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(self.to_json(), encoding="utf-8")


def _percentile(sorted_vals: list[float], p: float) -> float:
    if not sorted_vals:
        return 0.0
    k = (len(sorted_vals) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[f] * (c - k) + sorted_vals[c] * (k - f)


def _softclip_length(cigartuples) -> int:
    """Length of soft-clip portion of a CIGAR (SAM op-code 4)."""
    if not cigartuples:
        return 0
    total = 0
    for op, ln in cigartuples:
        if op == 4:
            total += ln
    return total


def compute_library_stats(
    bam_or_cram: str,
    reference_fasta: str | None = None,
    sample_size: int = 100_000,
    min_mapq_for_insert: int = 20,
) -> LibraryStats:
    """Sample up to ``sample_size`` reads from the BAM and infer library stats.

    Strategy: sample from a representative euchromatic window on chr1
    (30–60 Mb). Streaming from the *start* of chr1 only covers the first ~1 Mb —
    the subtelomeric repeat region where ~75% of reads are MAPQ 0 — which yields
    a false median MAPQ of 0 and biases the duplicate/insert stats. The mid-arm
    window is unique-mapping (median MAPQ ~60) and representative of the library.
    """
    open_kwargs = {}
    if bam_or_cram.endswith(".cram") and reference_fasta:
        open_kwargs["reference_filename"] = reference_fasta
    mode = "rc" if bam_or_cram.endswith(".cram") else "rb"
    sam = pysam.AlignmentFile(bam_or_cram, mode, **open_kwargs)

    insert_sizes: list[float] = []
    softclips: list[int] = []
    mapqs: list[int] = []
    n_total = 0
    n_supp = 0
    n_dup = 0

    # Choose a contig — prefer chr1, fall back to first reference.
    refs = list(sam.references)
    contig = "chr1" if "chr1" in refs else ("1" if "1" in refs else (refs[0] if refs else None))
    if contig is None:
        sam.close()
        return LibraryStats()

    # Representative euchromatic window on chr1 (avoids the subtelomeric repeats
    # at the contig start). Fall back to whole-contig for non-chr1 references.
    if contig in ("chr1", "1"):
        contig_len = sam.get_reference_length(contig)
        win = (30_000_000, 60_000_000) if contig_len and contig_len > 60_000_000 else (None, None)
        fetch_iter = sam.fetch(contig, win[0], win[1]) if win[0] is not None else sam.fetch(contig)
    else:
        fetch_iter = sam.fetch(contig)

    for read in fetch_iter:
        n_total += 1
        if n_total > sample_size * 4:    # safety cap on iteration
            break
        if read.is_supplementary or read.is_secondary:
            if read.is_supplementary:
                n_supp += 1
            continue
        if read.is_unmapped:
            continue
        if read.is_duplicate:
            n_dup += 1
            continue
        if read.mapping_quality >= min_mapq_for_insert:
            tlen = abs(read.template_length or 0)
            if 0 < tlen < 100_000:    # ignore degenerate / chimeric inferred sizes
                insert_sizes.append(float(tlen))
        softclips.append(_softclip_length(read.cigartuples or []))
        mapqs.append(int(read.mapping_quality))
        if len(insert_sizes) >= sample_size:
            break

    sam.close()

    if not insert_sizes:
        return LibraryStats(n_reads_sampled=len(mapqs),
                            pct_duplicate=n_dup / max(n_total, 1),
                            pct_supplementary=n_supp / max(n_total, 1))

    insert_sizes.sort()
    softclips.sort()
    mapqs.sort()

    median = statistics.median(insert_sizes)
    deviations = sorted(abs(x - median) for x in insert_sizes)
    mad = statistics.median(deviations)

    return LibraryStats(
        n_reads_sampled=len(insert_sizes),
        insert_median=median,
        insert_mad=mad,
        insert_p95=_percentile(insert_sizes, 0.95),
        insert_p99=_percentile(insert_sizes, 0.99),
        softclip_median=statistics.median(softclips) if softclips else 0.0,
        softclip_p95=_percentile(softclips, 0.95),
        mapq_median=statistics.median(mapqs) if mapqs else 0.0,
        mapq_p10=_percentile(mapqs, 0.10),
        pct_supplementary=n_supp / max(n_total, 1),
        pct_duplicate=n_dup / max(n_total, 1),
    )


def mapq_weight(mapq: int, full_weight_mapq: int = 30) -> float:
    """GRIDSS-style soft MAPQ weight.

    A read with MAPQ ≥ ``full_weight_mapq`` contributes 1.0. Below that, it
    contributes ``mapq / full_weight_mapq`` (linearly tapered). MAPQ = 0
    contributes 0.

    The strict-cutoff (``min_mapq``) is layered on top in the scanner — this
    weight only applies to reads that already passed the floor.
    """
    if mapq <= 0:
        return 0.0
    if mapq >= full_weight_mapq:
        return 1.0
    return mapq / float(full_weight_mapq)
