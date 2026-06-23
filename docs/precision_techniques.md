# Precision-improving techniques worth porting into quasarsv

Source-code survey of 7 competitor SV callers, mapped against quasarsv's
known false negatives:

* **NU-DHL-1** — documented BCL6-IGH (t(3;14)). Our pipeline emits zero
  evidence at any tier.
* **SU-DHL-9** — documented MYC-IGH (t(8;14)). Our pipeline emits zero
  evidence at any tier.
* **OCI-Ly19** — documented BCL2-IGH (t(14;18)). Our pipeline surfaces
  MYC-IGH only at T2 with PE=5 (truth pair absent).

Numbers below: F1 lift estimates are subjective until the benchmark
(currently running, see `output/benchmark/`) provides head-to-head
quasarsv-vs-each-caller comparisons against `cohort_truth.tsv`.

---

## Top 10 techniques to port (ranked by expected impact / port cost)

### 1. GRIDSS empirical-LLR variant scoring  *(high impact / medium cost)*

**Source:** `gridss/src/main/java/au/edu/wehi/idsv/model/EmpiricalLlrModel.java`
+ `MapqModel.java` (the fast variant)

**What it does:** Replaces hard tier thresholds (`SR>=3 AND PE>=2 AND ...`)
with a log-likelihood ratio scored against each library's empirical CIGAR
distribution and fragment-size distribution. Each piece of evidence
contributes a Phred-scaled LLR weighted by its MAPQ.

```
llr = log10( (prEgivenMR + prM * (prEgivenMV - prEgivenMR)) / prEgivenMR )
```

* `prEgivenMR` = probability of seeing this evidence under the reference
  (drawn from library empirical distribution).
* `prEgivenMV` = probability under the variant.
* `prM`       = MAPQ-weighted prior on the read mapping being correct.

**Why this lifts our F1:** OCI-Ly19's BCL2-IGH sits at PE=5 with no SR — our
fixed thresholds drop it. LLR scoring would assign it a real (non-zero)
score that surfaces it at T2/T1 if the library's discordant rate is low.
The whole tier system becomes adaptive to the sample's library.

**Port plan:**
1. Add `EmpiricalLibraryStats` (fragment-size histogram + CIGAR soft-clip
   length distribution) — computed once per BAM/CRAM at the start of
   `scan-cram`.
2. New `score_evidence_llr(evidence, stats)` helper in `merge.py`.
3. Replace the tier-rule cascade with a continuous score; tier thresholds
   become quantiles of the score distribution per sample.
4. Keep the rule-based path as a fallback when library stats can't be
   computed.

---

### 2. GRIDSS maximum-clique breakend interval clustering  *(high impact / medium cost)*

**Source:** `gridss/src/main/java/au/edu/wehi/idsv/model/Models.java`
+ `gridss/src/main/java/au/edu/wehi/idsv/graph/MaximumCliqueIntervalGraph.java`

**What it does:** Finds the maximum-weight interval that satisfies the
clique constraint: every pair of supporting evidence intervals must overlap.
Replaces our fixed 500 bp `pos // tolerance` buckets.

```java
Node fwd = maximalInterval(lgc, BreakendDirection.Forward,  bs, weights);
Node bwd = maximalInterval(lgc, BreakendDirection.Backward, bs, weights);
```

**Why this lifts our F1:** NU-DHL-1's BCL6-IGH evidence may be scattered
across multiple 500 bp buckets, none of them individually clearing the
`min_split_reads=3` threshold. A maximum-clique interval captures the
combined weight without artefacts of bucket boundaries.

**Port plan:**
1. New `merge/clique_cluster.py` implementing the interval-graph maximum-
   clique algorithm. (Pure-python prototype first; cython optimisation
   later if needed.)
2. Replace the `key = (sa_chrom, _bucket(pb, tol), ...)` clustering in
   `scanners/cram_scanner.py` with weighted-interval clustering.
3. Each cluster's `pos` becomes the weighted centroid of the clique
   (currently the median of pos_a_examples).

---

### 3. SvABA local-assembly-anchored breakpoints  *(high impact / high cost)*

**Source:** `svaba/src/svaba/AlignedContig.cpp` (esp. lines 60–280)
+ `svaba/src/svaba/BreakPoint.cpp`

**What it does:** Per locus, SvABA collects soft-clipped reads + their
mates, runs an SGA-style local assembler (`AssemblyEngine`), produces
contigs, then aligns each contig back to the reference with BWA-MEM.
Discontinuities in the contig alignment are breakpoints with single-base
resolution.

**Why this lifts our F1:** OCI-Ly19's PE=5 BCL2-IGH would, when locally
assembled, produce a contig spanning the junction with the IGH switch
sequence on one end and BCL2 on the other. The assembled contig adds a
new evidence type beyond PE+SR: `assembly_contigs`, which our schema
already has but currently always 0 for `forge_scan`.

**Port plan:**
1. Add `scanners/local_assembly.py` using `pyspoa` (Python POA assembler;
   no SGA dep) or shell out to `wtdbg2` for a heavier dep.
2. Trigger at any locus where SR + PE > threshold but neither alone is
   T1-eligible. Assemble + align contigs.
3. If a contig has a chimeric alignment, emit a new BreakpointCall with
   `assembly_contigs=1` + `precise=True`.

**Note:** Higher port cost because it adds a new external dep + an
assembly pipeline. Defer to v3 unless OCI-Ly19/NU-DHL-1 don't improve
from #1+#2 alone.

---

### 4. GRIDSS single-breakend support  *(medium impact / medium cost)*

**Source:** `gridss/src/main/java/au/edu/wehi/idsv/SingleReadEvidence.java`

**What it does:** When only one side of a putative SV has alignable reads
(e.g. the partner side is in a polyG attractor or unmapped), GRIDSS still
emits a "single breakend" call with the unaligned sequence in the ALT
field. Downstream filtering / annotation can then resolve the partner
later (e.g. via repeat-class lookup, BLAST, or longer reads).

**Why this lifts our F1:** Karpas-1106P's IGH side of t(3;14) is routed
to chr2:32916 polyG. A single-breakend call at BCL6 with the soft-clipped
sequence preserved would let downstream annotation match it against IG
switch consensus sequences.

**Port plan:**
1. New `BreakpointCall` variant: `chrom_b="."`, `pos_b=0`, but
   `soft_clip_consensus` populated.
2. Annotation step matches `soft_clip_consensus` against
   IG switch consensus sequences (S-µ, S-γ1-4, S-α1-2, S-ε) — these
   are ~30-mer GC-rich repeats with known motifs.
3. If the consensus matches a switch region, infer the IG side of the
   translocation directly.

---

### 5. Delly split-read realignment with BWA-MEM  *(medium impact / low cost)*

**Source:** `delly/src/junction.h` + `delly/src/cluster.h` (`SRAlignment`)

**What it does:** Delly takes soft-clipped reads, extracts the clipped
sequence, and realigns it with BWA-MEM Smith-Waterman. Confirms the SA-tag
partner with an independent aligner pass — catches cases where the SA
tag's aligner was wrong.

**Why this lifts our F1:** Our pysam scanner trusts the BAM's SA tag.
Some BAMs from older pipelines have stale or missing SA tags. A
ground-truth realignment of soft-clipped sequence catches those.

**Port plan:**
1. `scanners/sa_realigner.py` — given a soft-clipped read, run BWA-MEM
   on the clipped sequence (subprocess to the `base-bio` conda env).
2. Compare the realigner's hit against the SA tag. If they disagree
   or the SA tag is missing, use the realigner's hit.
3. Add this as a configurable step in `cram_scanner.scan_cram`.

---

### 6. Adaptive fragment-size threshold  *(medium impact / low cost)*

**Source:** Delly `cluster.h::libraryStat` + GRIDSS `IdsvSamFileMetrics`

**What it does:** Per BAM, compute the actual insert-size distribution
(median + MAD or kde). A pair is "discordant" if its insert size is
> 5× MAD from the median (Delly default), not a fixed 10 kb threshold.

**Why this lifts our F1:** Our `discordant_min_distance=10000` may be too
strict for some libraries and too lax for others. PCR-free libraries with
600 bp fragments would treat ~3 kb pairs as concordant under a 10 kb
threshold, missing nearby SVs.

**Port plan:**
1. New `scanners/library_stats.py` — sample 100k reads per BAM, compute
   median + MAD of insert size + soft-clip-length distribution.
2. Make `ScannerConfig.discordant_min_distance` automatically derived
   from `(median + 5 * MAD)` if not user-set.
3. Cache per-BAM stats in `output/<sample>/library_stats.json`.

---

### 7. TIDDIT DBSCAN clustering  *(medium impact / low cost)*

**Source:** `TIDDIT/tiddit/DBSCAN.py` + `tiddit/tiddit_cluster.pyx`

**What it does:** Uses DBSCAN (density-based clustering) on the read
endpoints instead of fixed-bucket clustering. Two reads cluster together
if they're within `eps` of each other AND there are at least `min_samples`
within `eps` of one of them. Handles arbitrary density without grid
artefacts.

**Why this lifts our F1:** Similar to GRIDSS max-clique but cheaper to
implement (sklearn has DBSCAN out of the box). The IG switch regions are
dense with breakpoints — DBSCAN would isolate the true IG-driver junction
from background switch-region rearrangement noise.

**Port plan:**
1. `from sklearn.cluster import DBSCAN` in `merge.py`.
2. For each `(chrom_a, chrom_b)` chromosome pair, run DBSCAN on the
   joint (pos_a, pos_b) coordinates.
3. eps = library median fragment size; min_samples = 2.

---

### 8. Lumpy probability-vector breakpoint  *(medium impact / medium cost)*

**Source:** `lumpy-sv/src/lumpy/SV_BreakPoint.h` + `SV_Pair.h`

**What it does:** Each breakpoint is represented as a probability vector
over an interval — not a single point. Evidence from PE / SR / soft-clip
each contributes its own probability distribution; combination is the
element-wise product (in log-space).

**Why this lifts our F1:** The maximum-likelihood breakpoint position
from the combined distribution is more accurate than median/centroid.
Useful when reporting precise breakpoints to clinicians (matters for
primer design in confirmatory PCR).

**Port plan:** Defer — orthogonal to recall improvement. Useful later
for reporting accuracy.

---

### 9. Manta candidate-junction graph traversal  *(high impact / high cost)*

**Source:** `manta/src/c++/lib/applications/EstimateSVLoci/` (graph
construction) + `manta/src/c++/lib/manta/SVLocusGraph.cpp` (traversal)

**What it does:** Builds a graph where nodes = genomic intervals, edges =
weighted by PE+SR evidence. Connected components correspond to multi-
breakpoint events (e.g. chained translocations, complex inversions).

**Why this lifts our F1:** NU-DHL-1's BCL6-IGH might be a complex
rearrangement (cell-line drift can introduce additional breaks). A
graph-traversal step captures the multi-hop event that our pairwise
clustering can't see.

**Port plan:**
1. Build a `networkx.Graph` post-merge: nodes = `(chrom, pos // 1kb)`,
   edges = `(node_a, node_b, weight=evidence)`.
2. Find connected components — each gets a "component_id" tagged onto
   every FusionCall.
3. Use component-level aggregation: a multi-edge component with low
   per-edge but high cumulative evidence gets a tier boost.

**Note:** High port cost; major architectural shift. Defer to v3.

---

### 10. MAPQ-as-weight instead of MAPQ-filter  *(low impact / low cost)*

**Source:** GRIDSS `MapqModel.java` (the entire 50-line file)

**What it does:** Where we hard-filter `read.mapping_quality < 20`,
GRIDSS uses MAPQ as a continuous weight: `phredOr(mapq1, mapq2)`.
A pair of mapq=15 reads contributes more than zero.

**Why this lifts our F1:** Reads near polyG attractors often have
mapq=0-10 due to multi-mapping. Hard-filtering loses 100% of their signal;
weighting still captures ~30% of it. Could surface NU-DHL-1's lost reads.

**Port plan:**
1. `cram_scanner.scan_cram`: change `if read.mapping_quality < min_mapq:
   continue` to keep low-mapq reads but record `mapq_weight = phred_to_pr(mapq)`.
2. Cluster SR/PE counts become weighted sums: `sum(mapq_weight)` not
   `len(reads)`.
3. Tier thresholds become weighted-sum thresholds.

---

## Per-tool summaries

**Manta** (`tools_src/manta/`) — Illumina's flagship. Graph-
based candidate generation, then per-candidate scoring with somatic and
diploid models. Strengths: assembly-anchored breakpoint refinement,
chromosome-depth-aware filtering. Weaknesses: slow on PMBL samples
(~30-90 min targeted on a 60 GB CRAM). Best ideas to port: graph
traversal (#9), assembly anchoring (subset of #3).

**GRIDSS** (`tools_src/gridss/`) — most algorithmically
mature. Single-breakend support, empirical-LLR scoring, maximum-clique
breakend intervals. Java-based, very slow (~6 hr per WGS). Best ideas:
techniques #1, #2, #4, #10. The single biggest leverage point in this
survey.

**Delly** (`tools_src/delly/`) — efficient C++. PE clustering
with read-depth integration. Strengths: speed, library-statistics-aware
thresholds. Weaknesses: less sensitive to SR-only events. Best ideas:
#5 (BWA-MEM realignment), #6 (library stats).

**SvABA** (`tools_src/svaba/`) — local-assembly-first design.
Slowest of the bunch (~3-4 hr per WGS) but produces single-base-resolution
breakpoints with assembled contigs. Best idea: #3 (local assembly).

**TIDDIT** (`tools_src/TIDDIT/`) — Python+Cython. DBSCAN
clustering, coverage-deviation-driven calling. Lightweight. Best idea:
#7 (DBSCAN).

**Lumpy** (`tools_src/lumpy-sv/`) — probability-vector
breakpoints. First mover in evidence-likelihood integration but largely
superseded by GRIDSS. Best idea: #8 (probability vectors) for breakpoint
reporting accuracy.

**SViCT** (`tools_src/svict/`) — cancer-fusion-specialist.
Low-VAF SV calling, soft-clip clustering for fusion breakpoints. Sparse
source repo (124 KB); main techniques are restatements of #5 and #7
specialised for fusion contexts. No unique technique worth porting.

---

## Coverage matrix — which techniques address which FN

| Technique | NU-DHL-1 (BCL6-IGH) | SU-DHL-9 (MYC-IGH) | OCI-Ly19 (BCL2-IGH) | General |
|---|---|---|---|---|
| #1 Empirical LLR     | ✓ | ✓ | **✓ (likely fix)** | ✓ |
| #2 Max-clique cluster| **✓** | **✓** | ✓ | ✓ |
| #3 Local assembly    | ✓ | ✓ | **✓** | ✓ |
| #4 Single-breakend   | ✓ |   |   | ✓ (PMBL polyG) |
| #5 BWA-MEM realign   | ✓ | ✓ |   | ✓ |
| #6 Adaptive insert   |   |   | ✓ | ✓ |
| #7 DBSCAN            | ✓ | ✓ |   | ✓ |
| #8 Prob vectors      |   |   |   | (precision only) |
| #9 Graph traversal   | **✓** | ✓ |   | ✓ (complex SVs) |
| #10 MAPQ weighting   | ✓ | ✓ |   | ✓ |

Bold = the technique most likely to be the difference-maker for that FN.

---

## Recommended port order (sequence + measure)

1. **#10 MAPQ-as-weight** — single-file change, 1 hour. Re-run benchmark.
2. **#6 Adaptive insert threshold** — add `library_stats.py`, 2 hours.
   Re-run benchmark.
3. **#7 DBSCAN clustering** — swap `merge.py` clustering, 2 hours.
   Re-run benchmark.
4. **#1 Empirical-LLR scoring** — new `EmpiricalLibraryStats` infrastructure,
   half day. Replace tier rules. **The big one.**
5. **#5 BWA-MEM realignment** — new scanner, half day. Optional unless
   the above leave residual FNs.
6. **#3 Local assembly** — only if FNs remain after the above.
7. **#2 Max-clique clustering** — only if DBSCAN (#7) doesn't capture the
   non-uniform-density cases.

Stop after every step, run `quasarsv benchmark` against
`cohort_truth.tsv`, and ship the change only if F1 moves forward without
regressing existing wins.

---

## Open question — read-data discriminator for IG-partner in PMBL

The PMBL precision regression (Karpas-1106P, U2940) showed that no
read-data signal cleanly favours IGH over IGL/IGK. The competitors don't
have a magic answer either — Manta and GRIDSS at chr2:32916 would
similarly see a flat SA-tag distribution. The actual fix space:

1. Match soft-clipped sequence against IG switch consensus motifs
   (#4 single-breakend + switch consensus lookup).
2. Use a cytogenetic prior (already done — B-cell lineage filter).
3. Long-read confirmation (out of scope until ONT/PacBio coverage exists).

None of these come for free from the surveyed callers.
