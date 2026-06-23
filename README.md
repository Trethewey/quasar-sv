<p align="center">
  <img src="brand/png/quasar-banner-dark.png" alt="Quasar" width="720">
</p>

A structural-variant caller built for lymphoma.

Quasar reads aligned sequencing data (BAM or CRAM) and finds the
translocations that matter clinically in lymphoma — the canonical IGH
partner fusions, BCL6 rearrangements, NPM1-ALK, and the rest. It writes
the calls to a VCF and a per-sample clinical report.

## Why

The general-purpose SV callers — Manta, GRIDSS, Delly, SvABA, TIDDIT —
are built for tumours in general, not lymphoma in particular. Three things
specifically:

1. **They miss IG-switch translocations.** Chimeric reads from t(3;14)
   BCL6-IGH and similar events route into a known GRCh38 reference
   artefact at chr2:32916 and never reach the true partner. Quasar's
   rescue layer reconstructs the underlying translocation from the
   artefact-mediated signal.

2. **They don't know which gene pairs are clinically meaningful.**
   Quasar carries a curated table of lymphoma canonical partners (IGH-BCL2,
   IGH-MYC, BCL6-IGH and the others) and tiers any annotated pair with
   non-trivial evidence to a clinical confidence level — even when a
   general-purpose caller would filter it as low-quality noise.

3. **They treat aSHM hotspot artefacts and V(D)J recombination as real
   SVs.** Quasar's classification layer recognises these patterns and
   demotes them out of the clinical tier unless multiple callers
   independently agree on them.

## Install

```bash
pip install quasar-sv
```

The CLI command is `quasar`. The Python package imports as `quasarsv`.

## Run

```bash
quasar call \
    --sample SAMPLE_ID \
    --bam SAMPLE.cram \
    --reference GRCh38.fasta \
    --output-dir output/SAMPLE_ID
```

Outputs in `output/SAMPLE_ID/`:

- **`SAMPLE_ID.fusions.tsv`** — every call with tier, partner genes,
  evidence summary, and QC flags
- **`SAMPLE_ID.fusions.vcf.gz`** — same calls as VCF 4.3
- **`brochure_SAMPLE_ID.html`** — per-sample clinical brochure

Tiers:

- **T1** — high confidence. A known canonical lymphoma partner pair
  with strong evidence, or a call backed by multiple independent signals.
- **T2** — moderate confidence. Annotated lymphoma partner with weaker
  evidence, or PASS in two or more callers.
- **T3** — surfaced for review; usually noise.

## Documentation

- `docs/algorithm_vignette.md` — full algorithm walkthrough
- `docs/benchmark_results.md` — head-to-head benchmark numbers + caveats
- `docs/panel_validation.md` — running on targeted panel BAMs
- `docs/quasar_vignette.docx` — printable summary

## Scope

Alpha release. Validated on a 14-sample lymphoma cell-line WGS cohort.
Held-out validation, broad-cohort generalisability, and patient-sample
performance are unmeasured. **Not for clinical decision-making without
site-specific validation.**

## Author

Chris Trethewey — [christrethewey.dev](https://christrethewey.dev/) — [github.com/Trethewey](https://github.com/Trethewey)

Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
