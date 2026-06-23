"""Read-level scanners — extract breakpoint evidence directly from BAM/CRAM.

Used when no upstream SV caller has been run. The scanner targets a configurable
set of loci (default: built-in lymphoma driver + IG loci) and identifies
candidate breakpoints from split reads (supplementary alignments, SA tag) and
discordant read pairs.

Requires pysam — install via the [bam] extra.
"""
from .cram_scanner import scan_cram, ScannerConfig, scan_to_breakpoint_calls
from .sa_aware import scan_artefacts_sa, SAScannerConfig
from .sa_chrom_inference import scan_artefacts_chrom_inference, ChromInferenceConfig

__all__ = ["scan_cram", "scan_to_breakpoint_calls", "ScannerConfig",
           "scan_artefacts_sa", "SAScannerConfig",
           "scan_artefacts_chrom_inference", "ChromInferenceConfig"]
