#!/usr/bin/env python3
"""One-shot cleanup: undo the LLR-promoted tier changes that were written to
disk while empirical-LLR was wired into apply_default_qc.

LLR promotion tagged each lifted call with `llr_promoted_<NEW>_from_<OLD>`.
We use that breadcrumb to restore the original tier and strip the flag.
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.model import read_fusion_calls_tsv, write_fusion_calls_tsv

ROOT = Path(__file__).resolve().parent.parent
RE_FLAG = re.compile(r"^llr_promoted_(T[123])_from_(T[123])$")


def main() -> int:
    total_calls = 0
    total_reverted = 0
    for tsv in sorted(glob.glob(str(ROOT / "output/wgs_cohort/*/*.fusions.tsv"))):
        calls = read_fusion_calls_tsv(tsv)
        n_reverted = 0
        for c in calls:
            kept_flags = []
            original_tier = None
            for f in c.qc_flags:
                m = RE_FLAG.match(f)
                if m:
                    original_tier = m.group(2)
                else:
                    kept_flags.append(f)
            if original_tier is not None and original_tier != c.tier:
                c.tier = original_tier
                n_reverted += 1
            c.qc_flags = kept_flags
        total_calls += len(calls)
        total_reverted += n_reverted
        if n_reverted:
            write_fusion_calls_tsv(calls, tsv)
            print(f"  {Path(tsv).name}: reverted {n_reverted} of {len(calls)} calls")
    print(f"\nTotal reverted: {total_reverted} / {total_calls}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
