#!/usr/bin/env bash
# Sequentially scan every CRAM in the WGS cohort. One CRAM at a time
# (bandwidth-limited NAS-friendly). Skips samples whose fusions.tsv exists.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

PROJECT="$PROJECT_ROOT"
REF="$FF_REFERENCE"
OUT="$FF_OUTPUT_DIR/wgs_cohort"
mkdir -p "$OUT"
LOG=$OUT/_batch.log

cd "$PROJECT"

log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" | tee -a "$LOG"; }

# Build sample manifest from cohort_metadata.xlsx if not already present
MANIFEST=$OUT/_samples.tsv
if [[ ! -s "$MANIFEST" ]]; then
  log "Building sample manifest from cohort_metadata.xlsx"
  python3 scripts/build_wgs_samples_tsv.py \
    --metadata "$FF_COHORT_METADATA" \
    --cram-root "$FF_CRAM_ROOT" \
    --ref "$REF" --output "$MANIFEST"
fi
N_TOTAL=$(wc -l < "$MANIFEST")
log "$N_TOTAL samples in manifest"

idx=0
while IFS=$'\t' read -r sample cram ref; do
  idx=$((idx+1))
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue
  if [[ ! -f "$cram" ]]; then
    log "[$idx/$N_TOTAL] SKIP $sample — CRAM not found: $cram"
    continue
  fi
  if [[ ! -f "$ref" ]]; then
    log "[$idx/$N_TOTAL] SKIP $sample — reference not found: $ref"
    continue
  fi
  outdir=$OUT/$sample
  if [[ -f "$outdir/${sample}.fusions.tsv" ]]; then
    log "[$idx/$N_TOTAL] EXISTS $sample — already scanned, skipping"
    continue
  fi
  mkdir -p "$outdir"
  log "[$idx/$N_TOTAL] SCAN $sample"
  t0=$(date +%s)
  python3 -m quasarsv.cli scan-cram \
    --sample "$sample" \
    --bam "$cram" \
    --reference "$ref" \
    --output-dir "$outdir" \
    --min-split-reads 3 --min-discordant-pairs 5 \
    --pad-locus-bp 50000 \
    --skip-reports \
    >> "$outdir/_scan.log" 2>&1 || {
      log "[$idx/$N_TOTAL] FAIL $sample (exit $?) — see $outdir/_scan.log"
      continue
    }
  dt=$(( $(date +%s) - t0 ))
  log "[$idx/$N_TOTAL] DONE $sample in ${dt}s"
done < "$MANIFEST"

log "Batch scan complete. Building cohort report from per-sample TSVs."

# Stitch every per-sample TSV into a single cohort table, then render dashboards.
python3 - <<EOF
import glob, sys
sys.path.insert(0, "src")
from quasarsv.model import read_fusion_calls_tsv, write_fusion_calls_tsv
all_calls = []
for p in sorted(glob.glob("$OUT/*/*.fusions.tsv")):
    all_calls.extend(read_fusion_calls_tsv(p))
print(f"merged {len(all_calls)} calls from per-sample TSVs")
write_fusion_calls_tsv(all_calls, "$OUT/cohort.fusions.tsv")
from quasarsv.reports import (
    write_brochure, write_cohort_dashboard, write_validation_report,
)
samples = sorted({c.sample for c in all_calls})
from quasarsv.metadata import load_cohort_metadata_xlsx, build_metadata_index
import os
COHORT_META = os.environ.get("FF_COHORT_METADATA", "")
try:
    meta_index = build_metadata_index(load_cohort_metadata_xlsx(COHORT_META)) if COHORT_META else {}
except Exception:
    meta_index = {}
for s in samples:
    sc = [c for c in all_calls if c.sample == s]
    write_brochure(s, sc, f"$OUT/{s}/brochure_{s}.html",
                   metadata=meta_index.get(s))
write_cohort_dashboard(all_calls, "$OUT/cohort_dashboard.html",
                       metadata_xlsx=COHORT_META if COHORT_META else None)
write_validation_report(all_calls, "$OUT/validation_report.html")

# Flat tabular cohort summary, one row per sample (for Excel review)
from collections import Counter
import csv
rows = []
for s in samples:
    sc = [c for c in all_calls if c.sample == s]
    tc = Counter(c.tier for c in sc)
    known_t1 = sorted({f"{c.gene_a}-{c.gene_b}" for c in sc
                        if c.tier == "T1" and c.known_partner})
    inferred = sorted({f"{c.gene_a}-{c.gene_b}" for c in sc
                       if "inferred_via_artefact_rescue" in c.qc_flags
                       and c.tier in ("T1", "T2")})
    drivers_t12 = sorted({c.driver_locus for c in sc
                          if c.driver_locus and c.tier in ("T1", "T2")})
    m = meta_index.get(s)
    rows.append({
        "sample": s,
        "cell_line": m.cell_line if m else "",
        "cohort": m.cohort if m else "",
        "ashm_expected": m.ashm_expected if m else "",
        "coverage_x": f"{m.coverage:.1f}" if m and m.coverage else "",
        "project": m.project if m else "",
        "T1": tc.get("T1", 0),
        "T2": tc.get("T2", 0),
        "T3": tc.get("T3", 0),
        "known_t1_partners": "; ".join(known_t1),
        "inferred_ig_partner_pairs": "; ".join(inferred),
        "drivers_with_T1_T2_hits": "; ".join(drivers_t12),
    })
with open("$OUT/cohort_summary.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"cohort_summary.tsv written: {len(rows)} rows")
print("rendered cohort dashboard + per-sample brochures + validation report")
EOF

log "Reports written under $OUT/"
