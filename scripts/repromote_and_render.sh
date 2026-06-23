#!/usr/bin/env bash
# Re-apply QC + canonical-partner promotion to every per-sample TSV and
# re-render brochures, cohort dashboard, validation report, and flat
# cohort summary. Does not re-scan any CRAM.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

PROJECT="$PROJECT_ROOT"
OUT="$FF_OUTPUT_DIR/wgs_cohort"
cd "$PROJECT"

python3 - <<EOF
import glob, sys
sys.path.insert(0, "src")
from quasarsv.model import read_fusion_calls_tsv, write_fusion_calls_tsv
from quasarsv.qc import apply_default_qc
from quasarsv.annotate import annotate_calls
from quasarsv.metadata import (
    load_cohort_metadata_xlsx, build_metadata_index, lineage_index_from_metadata,
)
from quasarsv.reports import write_brochure, write_cohort_dashboard, write_validation_report
from collections import Counter
import csv

try:
    meta_items = load_cohort_metadata_xlsx("${FF_COHORT_METADATA}")
    meta_index = build_metadata_index(meta_items)
    lineage_index = lineage_index_from_metadata(meta_items)
except Exception:
    meta_index = {}
    lineage_index = {}

all_calls = []
for p in sorted(glob.glob("$OUT/*/*.fusions.tsv")):
    calls = read_fusion_calls_tsv(p)
    # rerun annotation + qc (idempotent — qc flags dedupe; rescue handles existing)
    annotate_calls(calls)
    # Remove prior synthetic rescue entries before re-running (they're identifiable by qc_flag)
    calls = [c for c in calls if "inferred_via_artefact_rescue" not in c.qc_flags
             or c.callers_supporting != ["forge_scan_rescue"]]
    # Strip flags so qc can rebuild them
    for c in calls:
        c.qc_flags = [f for f in c.qc_flags if f not in
                      ("builtin_artefact_locus", "recurrent_artefact", "short_range",
                       "promoted_known_partner", "ig_partner_ambiguous")]
    apply_default_qc(calls, sample_lineage=lineage_index)
    write_fusion_calls_tsv(calls, p)
    all_calls.extend(calls)
print(f"re-tiered {len(all_calls)} calls across per-sample TSVs")

# Re-render per-sample brochures
samples = sorted({c.sample for c in all_calls})
for s in samples:
    sc = [c for c in all_calls if c.sample == s]
    write_brochure(s, sc, f"$OUT/{s}/brochure_{s}.html",
                   metadata=meta_index.get(s))
write_cohort_dashboard(all_calls, "$OUT/cohort_dashboard.html",
                       metadata_xlsx="${FF_COHORT_METADATA}")
write_validation_report(all_calls, "$OUT/validation_report.html")
write_fusion_calls_tsv(all_calls, "$OUT/cohort.fusions.tsv")

# Refreshed flat summary with canonical_partners_any_tier column
rows = []
for s in samples:
    sc = [c for c in all_calls if c.sample == s]
    tc = Counter(c.tier for c in sc)
    canon = sorted({f"{c.gene_a}-{c.gene_b}|{c.tier}" for c in sc
                    if c.event_class == "IG_driver_canonical"
                    and c.gene_a and c.gene_b})
    somatic_t12 = sum(1 for c in sc
                     if c.event_class in ("IG_driver_canonical", "driver_driver")
                     and c.tier in ("T1", "T2"))
    novel_t12 = sum(1 for c in sc
                    if c.event_class in ("IG_driver_novel", "driver_intergenic")
                    and c.tier in ("T1", "T2"))
    physio_n = sum(1 for c in sc
                   if c.event_class in ("IG_intra", "IG_IG"))
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
        "somatic_T1_T2": somatic_t12,
        "novel_T1_T2": novel_t12,
        "physiological_IG_TR": physio_n,
        "canonical_partners_any_tier": "; ".join(canon),
        "drivers_with_T1_T2_hits": "; ".join(drivers_t12),
    })
with open("$OUT/cohort_summary.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"cohort_summary.tsv refreshed ({len(rows)} rows)")
EOF
echo "done"
