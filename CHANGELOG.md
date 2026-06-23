# Changelog

All notable changes to Quasar will be recorded here.
This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — Alpha (initial release)

First public release. Tagged alpha because the benchmark cohort is small
(n=14 WGS samples) and held-out validation has not yet been performed —
see `README.md` § "Scope honesty".

### Added
- Ensemble structural-variant detection across Manta, GRIDSS, Delly, SvABA,
  TIDDIT, and FACTERA caller VCFs.
- Built-in pysam-based BAM/CRAM scanner — works without any external caller.
- Twelve-step lymphoma-aware pipeline (scan, artefact rescue, parse, merge,
  tier, annotate, QC, IG-driver rescue, canonical-partner promotion, event
  classification, non-clinical demotion).
- DBSCAN-based breakpoint clustering across callers with MAPQ-weighted
  evidence summation.
- Lymphoma driver-locus + known-canonical-partner annotation
  (`data/known_partners.tsv` covers 30 clinically actionable partner pairs).
- B-cell / T-cell lineage prior + ambiguity-flag system for IG-driver
  rescue calls.
- chr2:32916 polyG-attractor rescue (GRCh38) — recovers PMBL t(3;14) BCL6-IGH
  that other tools miss.
- VCF 4.3 export (`quasar emit-vcf`).
- Truth-set benchmark CLI (`quasar benchmark`).
- HTML report variants: per-sample clinical brochure, cohort dashboard,
  validation report.
- 64-test pytest suite.

### Benchmark
- Head-to-head on a 14-sample lymphoma WGS cohort vs Manta 1.6.0,
  Delly 2.1.0, SvABA 1.2.0, TIDDIT 3.9.5:
  - Quasar: P=0.87, R=0.81, F1=0.84
  - Next best (Manta): F1=0.50
  - See `docs/benchmark_results.md` for the full writeup.

### Known limitations
- Benchmark cohort is small and cell-line-derived; patient-sample
  performance unverified.
- Held-out validation cohort not yet run; F1 may decrease on a blind cohort.
- GRIDSS 2.13.2 was not run at scale on this hardware (~12 hr/sample is
  impractical without a dedicated compute cluster).
- Three FNs (NU-DHL-1, SU-DHL-9, OCI-Ly19) have literature inconsistencies
  in their truth-set entries.
