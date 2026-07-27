# Changelog

All notable changes to Quasar will be recorded here.
This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — Beta (initial release)

First public release. Tagged beta because the benchmark cohort is small
(n=14 WGS samples) and held-out validation has not yet been performed —
see `README.md` § "Scope and limitations".

### Added
- Ensemble structural-variant detection across Manta, GRIDSS, Delly, SvABA,
  TIDDIT, and FACTERA caller VCFs.
- Built-in pysam-based BAM/CRAM scanner — works without any external caller.
- Lymphoma-aware pipeline (scan, parse, merge, tier, annotate, QC,
  artefact-masked-breakend annotation, canonical-partner promotion, event
  classification, non-clinical demotion).
- DBSCAN-based breakpoint clustering across callers with MAPQ-weighted
  evidence summation.
- Lymphoma driver-locus + known-canonical-partner annotation
  (`data/known_partners.tsv` covers 30 clinically actionable partner pairs).
- chr2:32,916 poly-G attractor handling (GRCh38): adapter/poly-G clips are
  filtered at source, and a call whose partner breakend falls in a masked
  artefact region is marked `partner_undetermined` and capped at a review
  tier rather than being assigned a partner.
- Interchromosomal calls without discordant support are demoted: a real
  translocation puts read pairs across the junction, whereas a repeat that
  cross-maps a clipped read produces split reads and no pair.
- VCF 4.3 export (`quasar emit-vcf`).
- Truth-set benchmark CLI (`quasar benchmark`).
- HTML report variants: per-sample clinical brochure, cohort dashboard,
  validation report.
- 87-test pytest suite.

### Validation
- Scored on a 14-sample lymphoma cell-line WGS cohort against a truth set
  corrected against primary literature, and against an independent read-level
  junction oracle built directly from the CRAMs.
- Clinically actionable tier (T1): precision 1.000, recall 0.667, zero false
  positives — including zero on the two confirmed negative controls.
- Every true positive is backed by a direct read-level junction; none is
  attributable to a canonical-partner lookup.
- No head-to-head comparison is published. See "Known limitations".

### Known limitations
- Benchmark cohort is small and cell-line-derived; patient-sample
  performance unverified.
- Held-out validation cohort not yet run; performance may drop on a blind cohort.
- **No head-to-head comparison against other callers exists.** An earlier one
  was withdrawn: every external caller in it had been given region-restricted
  input by this project's own harness (Manta via `--callRegions`, SvABA via a
  pre-subset targeted BAM), so none could call an inter-chromosomal event, and
  the reported scores measured our configuration rather than their performance.
  Fair genome-wide re-runs are in progress. No comparative claim will be made
  until they finish and every caller is scored under identical rules.
- GRIDSS 2.13.2 has not been scored. It is slow rather than unusable
  (assembly + BWA realignment of every soft clip; hours per WGS sample), and
  is queued behind the other callers.
- Recall at T1 is 0.667: four confirmed events surface at T2 rather than T1,
  all with discordant support in the 5-9 range. Tier thresholds have not been
  tuned to the truth set, deliberately — tuning to it is what made the previous
  benchmark circular.
- `in_frame` is never populated: the tool has no transcript model. It requires
  a GTF that is not bundled.
- Four samples are excluded from scoring as unscoreable (contested identity or
  no named truth pair) and are reported as such rather than silently dropped.
