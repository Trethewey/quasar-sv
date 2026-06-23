"""Print a one-line-per-sample summary of all per-sample TSVs currently on disk.

Usage:
    python3 scripts/cohort_summary_so_far.py [--cohort-dir output/wgs_cohort]
                                              [--metadata <cohort_metadata.xlsx>]
"""
from __future__ import annotations

import argparse
import glob
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.model import read_fusion_calls_tsv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-dir", default="output/wgs_cohort")
    import os
    ap.add_argument("--metadata", default=os.environ.get("FF_COHORT_METADATA", ""))
    args = ap.parse_args()

    try:
        from quasarsv.metadata import (
            load_cohort_metadata_xlsx, build_metadata_index,
        )
        meta_index = build_metadata_index(load_cohort_metadata_xlsx(args.metadata))
    except Exception as e:
        print(f"# metadata load failed ({e}); proceeding without", file=sys.stderr)
        meta_index = {}

    header = f"{'sample':<32} {'cell_line':<16} {'cohort':<12} {'aSHM':<5} "\
             f"{'T1':>4} {'T2':>5} {'known partners':<32} {'putative IG-driver'}"
    print(header)
    print("-" * len(header))

    pattern = f"{args.cohort_dir}/*/*.fusions.tsv"
    paths = sorted(glob.glob(pattern))
    for p in paths:
        sname = Path(p).parent.name
        calls = read_fusion_calls_tsv(p)
        tc = Counter(c.tier for c in calls)
        known = sorted({f"{c.gene_a}-{c.gene_b}" for c in calls
                        if c.tier == "T1" and c.known_partner and c.gene_a and c.gene_b})
        inferred = sorted({f"{c.gene_a}-{c.gene_b}" for c in calls
                           if "inferred_via_artefact_rescue" in c.qc_flags
                           and c.tier in ("T1", "T2") and c.gene_a and c.gene_b})
        m = meta_index.get(sname)
        cell = (m.cell_line if m else "")[:15]
        cohort = (m.cohort if m else "")[:11]
        ashm = (m.ashm_expected if m else "")[:4]
        print(f"{sname:<32} {cell:<16} {cohort:<12} {ashm:<5} "
              f"{tc.get('T1', 0):>4} {tc.get('T2', 0):>5} "
              f"{', '.join(known)[:31]:<32} {', '.join(inferred)[:40]}")
    print(f"\n{len(paths)}/26 samples completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
