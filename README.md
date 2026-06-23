# Quasar

Lymphoma structural-variant and fusion caller. Reads aligned reads
(BAM/CRAM) and writes VCF.

Built for lymphoma diagnostics: detects the canonical translocations
(IGH-BCL2, IGH-MYC, BCL6-IGH and the rest), handles IG-switch-region
artefacts in GRCh38, and tiers calls by clinical confidence.

Status: alpha (v0.1.0).

## Install

```bash
pip install quasar-sv
```

CLI command is `quasar`. Python package is `quasarsv`.

## Use

```bash
quasar scan-cram \
    --sample SAMPLE_ID \
    --bam aligned.cram \
    --reference GRCh38.fasta \
    --output-dir output/SAMPLE_ID
```

Outputs in `output/SAMPLE_ID/`:

- `SAMPLE_ID.fusions.tsv` — schema-of-record (one row per call, tiered T1/T2/T3)
- `SAMPLE_ID.fusions.vcf.gz` — VCF 4.3 export
- `brochure_SAMPLE_ID.html` — per-sample clinical report

## What it does

Twelve-step pipeline. Full walkthrough in `docs/algorithm_vignette.md`.

- pysam scan of driver loci, IG/TR loci, and known reference artefact
  hotspots — split-read and discordant-pair evidence
- DBSCAN clustering of breakpoint pairs
- Tiering on independent evidence types (T1 / T2 / T3)
- Lymphoma annotation: driver-locus tagging + canonical-partner
  promotion from `data/known_partners.tsv`
- chr2:32916 polyG-attractor rescue — recovers IG-driver translocations
  whose chimeric reads route into the GRCh38 reference artefact
- Demotion of intra-gene single-caller T1 calls (V(D)J / aSHM-hotspot
  noise) unless multi-caller PASS preserves them

## Scope (alpha)

Tested on a 14-sample lymphoma WGS cell-line cohort.
See `docs/benchmark_results.md` for the benchmark numbers and the
full list of caveats (held-out validation not yet run, patient-sample
performance unmeasured, GRIDSS not run at scale).

Not for clinical decision-making without site-specific validation.

## Author + license

Chris Trethewey · [christrethewey.dev](https://christrethewey.dev/) ·
[github.com/Trethewey](https://github.com/Trethewey)

Apache-2.0. See `LICENSE` and `NOTICE`.
