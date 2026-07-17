# quasarsv — algorithm vignette

This is the working specification of what quasarsv actually does, end to
end, for a single sample. The format is sequential: every step lists the
inputs it consumes, the operation it performs, and the outputs it produces.
Source-file references are the authoritative implementation.

```
            ┌──────────────────────────────────────────────────┐
   inputs   │  BAM / CRAM (any of)                             │
            │  Manta / GRIDSS / Delly / SvABA / TIDDIT VCFs    │
            │  FACTERA panel TSV                               │
            └──────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   1  │  Per-locus pysam scan                                 │
      │   • driver + IG/TR + artefact loci (~40 windows)      │
      │   • SR (split-read via SA tag) + PE (discordant pair) │
      │   • emits BreakpointCall, caller="quasar"         │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   2  │  SA-tag scan of artefact loci                         │
      │   • clusters SA-target positions at chr2:32916        │
      │   • emits BreakpointCall, caller="quasar_sa"      │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   3  │  Chromosome-level SA inference (opt-in)               │
      │   • aggregates SA-target chromosomes at the artefact  │
      │   • caller="quasar_chrom_sa"                      │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   4  │  External VCF parsing (when present)                  │
      │   parsers/{manta,gridss,delly,svaba,tiddit,factera}.py│
      │   normalise to BreakpointCall                         │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   5  │  Evidence-level merger  (merge.py)                    │
      │   • cluster by (chrom_a, chrom_b, pos±250bp)          │
      │   • SUM independent evidence types, not caller votes  │
      │   • assign provisional tier (T1/T2/T3)                │
      │   • emits FusionCall                                  │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   6  │  Annotation  (annotate.py)                            │
      │   • gene_a / gene_b / region (exonic/up/downstream)   │
      │   • driver_locus label, known_partner flag            │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   7  │  QC flagging  (qc.py)                                 │
      │   • builtin artefact mask (chr2:32916 polyG)          │
      │   • recurrent-position flag (many partners, many chr) │
      │   • short-range intra-chr flag (aSHM-aware)           │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   8  │  Artefact-masked breakend annotation  (rescue.py)     │
      │   • partner_undetermined; no partner is invented      │
      │   • canonical-partner IG retention                    │
      │   • primary + ambiguous alternatives per top driver   │
      │   • fan-out control on weak shared signals            │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
   9  │  Known-canonical promotion  (promote.py)              │
      │   • known_partner pair + sufficient evidence → ≥T2    │
      │   • strong SR or PE → T1                              │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
  10  │  Event classification  (classify.py)                  │
      │   IG_intra / IG_IG / IG_driver_canonical /            │
      │   IG_driver_novel / driver_driver / driver_intra /    │
      │   driver_intergenic / intergenic                      │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
  11  │  Non-clinical T1 demotion  (classify.py)              │
      │   T1 → T3 for IG_intra/IG_IG/driver_intra/recurrent_  │
      │   artefact, unless multi-caller PASS or known partner │
      └───────────────────────────────────────────────────────┘
                              │
                              ▼
      ┌───────────────────────────────────────────────────────┐
  12  │  Outputs                                              │
      │   FusionCall TSV + JSON                               │
      │   per-sample HTML brochure                            │
      │   cohort dashboard                                    │
      │   validation report                                   │
      │   (optional) truth-set benchmark TSV                  │
      └───────────────────────────────────────────────────────┘
```

---

## 1 — Per-locus pysam scan
### Source: `src/quasarsv/scanners/cram_scanner.py`

For each of ~40 windows (driver loci from `data/lymphoma_loci.tsv` + IG/TR
loci + the chr2:32916 artefact), pysam-fetches every read overlapping the
window padded by 5 kb (`pad_locus_bp`).

For each read passing `MAPQ ≥ 20`, not duplicate / not secondary:

* **Split-read evidence** — if the read carries an `SA` tag, parse the
  supplementary alignment: `chrom`, `pos`, `strand`, `CIGAR`. The
  soft-clip side of the primary read determines the breakpoint orientation
  (leading clip → `-` strand, trailing clip → `+`).
* **Discordant-pair evidence** — if the read isn't soft-clipped but its
  mate maps to a different chromosome or `≥ 10 kb` away on the same
  chromosome, count it as PE.

Reads are clustered into per-(`mate_chrom`, `mate_pos // 500 bp`,
`strand_a`, `strand_b`) buckets. Each cluster gets a `BreakpointCall`
record with caller `quasar` and the cluster's representative position
(median of `pos_a_examples`), provided it meets
`min_split_reads ≥ 2 OR min_discordant_pairs ≥ 4`.

This is the only step that touches reads on big WGS. It runs over
random-access CRAM index (~30-90 s per CRAM for the lymphoma window set).

## 2 — SA-tag scan of artefact loci
### Source: `src/quasarsv/scanners/sa_aware.py`

Specifically targets the polyG attractor at chr2:32915800-32916800 (and any
other entries in `data/artefact_loci.tsv`).

Polymorphic-mapping artefacts absorb soft-clipped reads whose true partner
is somewhere else (typically an IG switch region). At the artefact, the
SA-tag points to the read's true other end. We:

1. Fetch every read overlapping the artefact (cap `max_reads_per_locus =
   200,000` since these loci attract millions).
2. For each read with an SA tag landing on a different chromosome, cluster
   by `(SA_chrom, SA_pos // 500, strand_a, strand_b)`.
3. Emit `BreakpointCall` records keyed on the TRUE partner (chrom_a)
   versus the artefact (chrom_b), caller `quasar_sa`.

This is the "reverse-direction" scanner. Note the geometry: at the attractor the
ALIGNED segment is the poly-G run and the soft-clip holds real sequence — the
mirror of the driver-side view — so the noise filter tests the aligned segment
here. Reads whose aligned segment is adapter/poly-G are rejected: their SA tag
names the junk read's origin locus, not a translocation partner. What survives
is a genuine chimera that happens to overlap the artefact window.

## 3 — Chromosome-level SA inference — REMOVED
### Source: `src/quasarsv/scanners/sa_chrom_inference.py`

Formerly opt-in via `scan-cram --chrom-sa-inference`: at each artefact locus it
tallied SA-tag chromosomes and recorded the dominant non-self chromosomes as
putative partners, emitting the **median of every SA position on that
chromosome** as the breakpoint coordinate and naming whichever gene lay nearest.

Removed. Two independent defects:

* The coordinate was invented. Those SA positions are unrelated to one another,
  so their median is an arbitrary point that need not lie near any junction —
  and a median that happens to land near IGH yields a confident-looking IGH
  partner.
* The input is not translocation signal. Reads inside the attractor are
  library-wide adapter/poly-G junk, so ranking their SA chromosomes ranks
  chromosomes roughly by size, and every large chromosome clears a 2% share.

The vignette's own note recorded the symptom — "the SA-tag distribution at
chr2:32916 is roughly proportional to baseline read depth and doesn't
preferentially point to chr14 (IGH)" — which is exactly what a signal-free
channel looks like. The function now raises `NotImplementedError`. Diffuse
partner signal is recovered from discordant mates, which carry a real
coordinate.

## 4 — External VCF parsing
### Source: `src/quasarsv/parsers/{manta,gridss,delly,svaba,tiddit,factera}.py`

Each parser normalises its caller's idioms (BND ALT syntax, FORMAT field
layout, INFO tags) into the unified `BreakpointCall`/`Evidence` schema.
Strand convention from VCF 4.3 BND: `t[p[` = forward-forward,
`t]p]` = forward-reverse, `[p[t` = reverse-forward, `]p]t` = reverse-
reverse.

Each parser emits a `caller` tag (`manta`, `gridss`, etc.) so the merger
can tell evidence types apart.

## 5 — Evidence-level merger
### Source: `src/quasarsv/merge.py` (`merge_caller_calls`, `TierThresholds`)

**This is the heart of the algorithm.** Naive caller-vote ensembling
double-counts: if Manta and GRIDSS both see the same split-read cluster,
that's one piece of evidence not two. quasarsv clusters across callers
by breakpoint position (±250 bp tolerance) and then sums *independent
evidence types* — split reads, discordant pairs, assembly contigs, soft
clips — and only THEN tiers on that.

For each merged cluster, a `FusionCall` is born with:

* `n_callers` = distinct callers contributing
* `callers_supporting` = list of caller names
* `n_evidence_types` = how many of {SR, PE, assembly, soft-clip} are present
* `split_reads`, `discordant_pairs`, `assembly_contigs`, `soft_clips` =
  max across contributing callers (avoids double-counting)
* `any_pass` = at least one caller's FILTER == PASS
* `precise` = any caller marked the breakpoint precise

Provisional tier (rewritten in step 8/11):

* **T1** = `(n_callers ≥ 2 AND any_pass AND split_reads ≥ 3 AND
            n_evidence_types ≥ 2 AND precise)`
          OR `(single caller with SR ≥ 10 AND PE ≥ 10 AND precise)`
          OR `(known canonical partner AND SR ≥ 5 OR PE ≥ 10)`
* **T2** = `(n_callers ≥ 2)` OR `(single PASS with SR ≥ 5 + assembly)`
          OR `(known canonical partner with SR + PE ≥ 5)`
* **T3** = everything else

The single-caller `quasar` path can reach T1 only via the SR ≥ 10 +
PE ≥ 10 rule or via the canonical-partner promotion in step 9. By design.

## 6 — Annotation
### Source: `src/quasarsv/annotate.py`

For each breakpoint side, looks up the gene at that position in
`data/lymphoma_loci.tsv` (with ±5 kb upstream/downstream pad). Sets:

* `gene_a`, `gene_b`
* `region_a`, `region_b` (one of `exonic_or_intronic`, `upstream`,
  `downstream`, `intergenic`)
* `driver_locus` = `"BCL6-IGH"` etc. when both sides are
  driver/IG-annotated
* `known_partner` = True if `(gene_a, gene_b)` is in
  `data/known_partners.tsv` (30 canonical lymphoma pairs); records the
  pair's `cytoband` and `disease` source

## 7 — QC flagging
### Source: `src/quasarsv/qc.py`

Three filters, advisory by default:

1. **`flag_builtin_artefact_loci`** — any breakpoint landing in
   `data/artefact_loci.tsv` is flagged `builtin_artefact_locus`. Tier is
   forced down to T3 (this is the only auto-downgrade — these are reliably
   not real). Step 8 marks these as partner_undetermined; nothing is
synthesised from them.
2. **`flag_recurrent_position_artefacts`** — a `(chrom, pos//500bp)`
   window with ≥10 distinct partners on ≥3 chromosomes is flagged
   `recurrent_artefact`. Not auto-downgraded; consumed by step 11.
3. **`flag_short_range_intrachr`** — intra-chr pairs `< 200 bp` apart
   are flagged `short_range`, UNLESS one side is in `ASHM_TARGETS`
   (`BCL6, BCL7A, BTG1, BTG2, MYC, PAX5, PIM1, CXCR4, IRF4, RHOH,
   ST6GAL1, SOCS1, REL, CIITA`) — clustered short-range breaks at aSHM
   targets are real B-cell biology.

## 8 — Artefact-masked breakend annotation
### Source: `src/quasarsv/rescue.py`

**This step used to synthesise driver-IG fusions. It no longer does, and the
reason matters more than the mechanism.**

The old logic: when step 7 flagged calls as `builtin_artefact_locus`, it assumed
the underlying biology was a real driver-IG translocation whose IG-switch side
had mismapped to the poly-G attractor, and reconstructed the pair — choosing the
IG partner from `known_partners.tsv` plus a B-cell lineage prior, with **no read
linking driver to partner**, then promoting canonical pairs to T1 with the
artefact-side read counts displayed as junction support.

That inference was invalid. Measured on the WGS validation cohort:

* **The channel carries no signal.** Every locus sheds reads to chr2:32,916 at
  ~200-280 split reads per 10k, rearranged or not. In Karpas-1106P, BCL6
  (204/10k, rearrangement claimed) sheds *less* than TP53 (240/10k) and NOTCH1
  (337/10k), neither of which is rearranged. "Both loci hit the same artefact"
  is true of every pair of loci in the genome.
* **The reads are not genomic.** The clipped segments are Illumina adapter
  read-through and 2-colour poly-G tails — identical at BCL6, IGH, TP53 and ATM.
  No threshold on sequencing chemistry can yield a partner.
* **The masking premise was false.** Discordant mates do not depend on junction
  sequence, so a poly-G attractor cannot hide them. MD903 shows a genuine
  BCL6-IGH at PE=16 by pairs alone. Karpas-1106P and U2940 — the two samples the
  rescue "uniquely caught" — have zero BCL6-IGH reads of either kind, and
  break-apart FISH (Dai 2015, PMID:26599546) shows BCL6 and IG-H/K/L germline in
  both. The translocation is absent, not masked.

What the step does now: mark any call whose partner breakend falls in a masked
artefact region with `partner_undetermined` and cap it at a review tier, so an
unresolvable breakend can never be presented as a resolved fusion. It invents
nothing and emits no synthetic calls. `rescue_ig_driver_pairs()` raises
`NotImplementedError`; the lineage prior is gone entirely, because its only
consumer was choosing a partner to name.

A real driver-IG translocation is recovered by the ordinary scanner path, from
reads that actually join the two loci.
## 9 — Known-canonical promotion
### Source: `src/quasarsv/promote.py`

Final pass per FusionCall: if `known_partner == True` and there's
non-trivial evidence:

* `SR ≥ 5` OR `PE ≥ 10` → tier becomes T1
* `SR + PE ≥ 5` → tier becomes T2

A call that survives step 7's mask and has a canonical gene pair WILL
reach at least T2. Tier never moves backward here — only up.

## 10 — Event classification
### Source: `src/quasarsv/classify.py`

Assigns `event_class` based on gene_a/gene_b roles in `lymphoma_loci.tsv`:

| Class | Both ends are... | Same gene? |
|---|---|---|
| `IG_intra`              | both IG/TR              | yes (V(D)J or switch) |
| `IG_IG`                 | both IG/TR              | no (usually noise)    |
| `IG_driver_canonical`   | one IG/TR + one driver  | pair in known_partners |
| `IG_driver_novel`       | one IG/TR + one driver  | pair NOT in known_partners |
| `IG_intergenic`         | one IG/TR + unannotated | — |
| `driver_driver`         | two drivers             | different genes |
| `driver_intra`          | two drivers             | same gene (aSHM duplicate) |
| `driver_intergenic`     | one driver + unannotated| — |
| `intergenic`            | neither annotated       | — |

`SOMATIC_CLINICAL = {IG_driver_canonical, driver_driver}` drives the
clinical KPI counts in the brochure / cohort dashboard.

## 11 — Non-clinical T1 demotion
### Source: `src/quasarsv/classify.py::demote_nonclinical_t1`

Single-caller `quasar` can hit T1 via the SR≥10 + PE≥10 path on
classes that aren't clinically actionable — V(D)J recombination at IG
loci (`IG_intra`), inter-IG noise (`IG_IG`), intra-driver aSHM
duplications (`driver_intra`), and `recurrent_artefact`-flagged windows.

This step demotes those T1 calls to T3 EXCEPT when:

* `known_partner == True` — canonical pair always preserved
* `n_callers ≥ 2 AND any_pass` — multi-caller agreement preserves T1 even
  for borderline classes

Adds `demoted_nonclinical_t1` to `qc_flags` when applied. Raises strict
precision by ~50 % on a typical sample without losing any canonical call.

## 12 — Outputs
### Source: `src/quasarsv/model.py` + `reports/`

* `<sample>.fusions.tsv` — the schema-of-record (32 columns, one row per
  FusionCall)
* `<sample>.fusions.json` — same payload for downstream tools
* `brochure_<sample>.html` — per-sample clinical brochure. Section order:
  sample provenance KPIs → detection summary → circos → canonical somatic
  translocations (per-class cards) → driver-driver fusions → putative
  novel → physiological IG/TR (top 15) → driver-locus hits → locus
  close-ups → QC summary.
* `cohort_dashboard.html` — across-sample view: per-translocation cards,
  per-sample KPI table, recurrent-rearrangement table.
* `validation_report.html` — per-caller Jaccard concordance, tier
  composition stack, replicate concordance (when supplied), PR curves
  (when `--truth-set` supplied).
* `cohort_summary.tsv` — flat one-row-per-sample table for Excel.

If invoked with `quasarsv benchmark`, additionally:

* `benchmark_<mode>.tsv` — per-sample matched / missed / FP-T1 + the
  aggregate TP/FP/FN + precision/recall/F1 line.

---

## Configuration surface (what the operator turns)

```python
ScannerConfig(min_mapq=20, discordant_min_distance=10_000,
              pos_tolerance=500, min_split_reads=2,
              min_discordant_pairs=4, pad_locus_bp=5_000)

MergeConfig(pos_tolerance=250, ...)

TierThresholds(...)   # the T1/T2/T3 rule constants

RescueConfig(
    unresolved_tier="T3",   # tier ceiling for an unresolvable breakend
)
# The partner-inference knobs are gone: they tuned a fabrication. No threshold
# on a signal-free channel can produce a valid partner assignment.
```

`--lineage`, `--metadata` and `--chrom-sa-inference` have been removed: the
first two fed only the partner inference, and the third fabricated a
breakpoint coordinate from a median of unrelated positions.

---

## Where it currently sits

Measured against `data/cohort_truth.tsv` — corrected against primary literature —
on the 14-sample cell-line WGS cohort, and against an independent read-level
junction oracle built directly from the CRAMs:

* **Clinically actionable (T1):** P = 1.000, R = 0.667, **F1 = 0.800**, with
  **zero false positives**, including zero on both confirmed negative controls.
* **Review list (T1+T2):** P = 0.375, R = 1.000, F1 = 0.545.
* **Detection vs lookup: 12 detected, 0 lookup-only.** Every true positive is
  backed by a read-level junction; none is attributable to a partner lookup.

The relaxed metric this section used to lead with is gone. `--relax-canonical-ig-
partner` let BCL6-IGL match a truth of BCL6-IGH, and a driver-only call match a
driver-IG truth; combined with a canonical-partner prior and a truth set encoding
the same canonical translocations, it scored the tool for identifying a driver
and looking its textbook partner up. It survives only as an explicitly labelled
side-metric and must never be the headline.

Four confirmed events surface at T2 rather than T1 (all with discordant support
of 5-9), which is the recall gap. Tier thresholds have deliberately NOT been
tuned to close it: tuning to the truth set is what made the previous benchmark
circular.

**No comparison against Manta / Delly / SvABA / GRIDSS / TIDDIT is published.**
An earlier one was withdrawn — every external caller had been run on
region-restricted input by this project's own harness, so their scores measured
our configuration, not their performance. Fair genome-wide re-runs are in
progress.
