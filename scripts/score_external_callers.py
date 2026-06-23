#!/usr/bin/env python3
"""Score external SV callers' VCFs against the cohort truth set.

For each (sample, tool) pair with a VCF in output/benchmark/<tool>/<sample>/:

  1. Parse the VCF through the matching quasarsv parser
     -> list[BreakpointCall]
  2. Wrap as a single-caller FusionCall list (no merge across callers — we are
     scoring each tool standalone). Each BreakpointCall is converted to a
     minimal FusionCall by the merge module's single-record path.
  3. Annotate gene_a/gene_b/known_partner via annotate_calls
  4. Apply the same default QC + tier logic quasarsv uses (so each tool is
     scored under matched conditions — promotion + classification active).
  5. Run quasarsv.benchmark.score_calls_against_truth against the packaged
     cohort_truth.tsv.

Outputs:
  output/benchmark/<tool>/<sample>/<sample>.<tool>.fusions.tsv
  output/benchmark/scores_<tool>.tsv         (per-tool per-sample TP/FP/FN)
  output/benchmark/scores_comparative.tsv    (quasarsv vs all callers)
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))

from quasarsv.annotate import annotate_calls
from quasarsv.benchmark import (
    builtin_truth_path, load_truth_set,
    score_calls_against_truth, write_benchmark_tsv,
)
from quasarsv.merge import MergeConfig, merge_caller_calls
from quasarsv.model import (
    FusionCall, BreakpointCall,
    write_fusion_calls_tsv,
)
from quasarsv.parsers import parse_any
from quasarsv.qc import apply_default_qc


TOOL_VCF = {
    "manta":  "results/variants/diploidSV.vcf.gz",
    "delly":  "{sample}.delly.vcf.gz",
    "svaba":  "{sample}.svaba.svaba.sv.vcf",
    "gridss": "{sample}.gridss.vcf.gz",
    "tiddit": "{sample}.tiddit.vcf",
}


def vcf_path_for(tool: str, sample: str, base: Path) -> Path | None:
    sub = base / tool / sample
    if not sub.exists():
        return None
    spec = TOOL_VCF.get(tool)
    if spec is None:
        return None
    p = sub / spec.format(sample=sample)
    return p if p.exists() else None


def score_one(tool: str, sample: str, base: Path,
              truth, sample_lineage: dict[str, str]):
    vcf = vcf_path_for(tool, sample, base)
    if vcf is None:
        return None

    # Parse VCF -> BreakpointCalls
    bps: list[BreakpointCall] = parse_any(str(vcf), tool, sample)
    if not bps:
        return {"vcf": str(vcf), "n_records": 0,
                "tp": 0, "fp": 0, "fn": 0, "skipped": True}

    # Merge into FusionCalls (single-caller "ensemble" of just this tool)
    per_caller = {tool: bps}
    calls = merge_caller_calls(per_caller, cfg=MergeConfig(pos_tolerance=250))
    annotate_calls(calls)
    apply_default_qc(calls, sample_lineage=sample_lineage)

    out_tsv = base / tool / sample / f"{sample}.{tool}.fusions.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_fusion_calls_tsv(calls, str(out_tsv))

    bench = score_calls_against_truth(
        calls, truth, tool_name=tool,
        relax_canonical_ig_partner=True,
        restrict_to_scored_samples=True,
    )
    return {
        "vcf": str(vcf),
        "n_records": len(bps),
        "n_fusion_calls": len(calls),
        "tp": bench.tp, "fp": bench.fp, "fn": bench.fn,
        "precision": bench.precision, "recall": bench.recall, "f1": bench.f1,
        "out_tsv": str(out_tsv),
        "bench": bench,
    }


def load_cohort_lineage_index():
    """Per-sample lineage from packaged cohort_truth.tsv cohort labels."""
    try:
        from quasarsv.metadata import (
            load_cohort_metadata_xlsx, lineage_index_from_metadata,
        )
        import os
        meta_path = os.environ.get("FF_COHORT_METADATA", "")
        if not meta_path:
            return {}
        items = load_cohort_metadata_xlsx(meta_path)
        return lineage_index_from_metadata(items)
    except Exception:
        return {}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=str(PROJECT / "output/benchmark"),
                   help="root containing per-tool VCFs")
    p.add_argument("--truth-set", default=str(builtin_truth_path()))
    p.add_argument("--tools", default="manta,delly,svaba,tiddit,gridss")
    p.add_argument("--samples",
                   help="optional comma-sep sample-id list; default = all in base/*/")
    args = p.parse_args()

    base = Path(args.base)
    truth = load_truth_set(args.truth_set)
    lineage = load_cohort_lineage_index()
    tools = [t.strip() for t in args.tools.split(",") if t.strip()]

    if args.samples:
        samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    else:
        samples_seen = set()
        for tool in tools:
            tdir = base / tool
            if tdir.exists():
                samples_seen.update(d.name for d in tdir.iterdir() if d.is_dir())
        # Always include samples whose truth entries are positive — gives a
        # meaningful quasarsv-only comparison even before tool VCFs exist.
        for t in truth:
            if t.is_positive:
                samples_seen.add(t.sample_id)
        samples = sorted(samples_seen)

    rows = []
    per_tool_aggregate = defaultdict(lambda: {"tp": 0, "fp": 0, "fn": 0,
                                              "n_samples": 0, "n_records": 0})
    for sample in samples:
        for tool in tools:
            r = score_one(tool, sample, base, truth, lineage)
            if r is None:
                continue
            rows.append({"sample": sample, "tool": tool, **{k: v for k, v in r.items()
                                                             if k != "bench"}})
            agg = per_tool_aggregate[tool]
            agg["tp"] += r["tp"]; agg["fp"] += r["fp"]; agg["fn"] += r["fn"]
            agg["n_samples"] += 1; agg["n_records"] += r.get("n_records", 0)

    # Add quasarsv's own per-sample scoring for comparison
    from quasarsv.model import read_fusion_calls_tsv
    ff_calls = []
    for sample in samples:
        tsv = PROJECT / "output/wgs_cohort" / sample / f"{sample}.fusions.tsv"
        if tsv.exists():
            ff_calls.extend(read_fusion_calls_tsv(str(tsv)))
    ff_bench = score_calls_against_truth(
        ff_calls, truth, tool_name="quasarsv",
        relax_canonical_ig_partner=True, restrict_to_scored_samples=True,
    )
    per_tool_aggregate["quasarsv"] = {
        "tp": ff_bench.tp, "fp": ff_bench.fp, "fn": ff_bench.fn,
        "n_samples": len(samples),
        "n_records": sum(1 for _ in ff_calls),
    }
    for sample, s in ff_bench.per_sample.items():
        rows.append({
            "sample": sample, "tool": "quasarsv",
            "vcf": "(internal)", "n_records": "-",
            "n_fusion_calls": len([c for c in ff_calls if c.sample == sample]),
            "tp": len(s.matched_pairs), "fp": s.fp_t1_count, "fn": len(s.missed_pairs),
            "precision": "", "recall": "", "f1": "",
            "out_tsv": str(PROJECT / "output/wgs_cohort" / sample
                           / f"{sample}.fusions.tsv"),
        })

    out_dir = base
    out_dir.mkdir(parents=True, exist_ok=True)
    long_path = out_dir / "scores_long.tsv"
    with open(long_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["sample", "tool", "n_records",
                                       "n_fusion_calls", "tp", "fp", "fn",
                                       "precision", "recall", "f1",
                                       "vcf", "out_tsv"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    comparative_path = out_dir / "scores_comparative.tsv"
    with open(comparative_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["tool", "n_samples", "n_records",
                    "tp", "fp", "fn", "precision", "recall", "f1"])
        order = ["quasarsv"] + [t for t in tools if t in per_tool_aggregate
                                    and t != "quasarsv"]
        for tool in order:
            agg = per_tool_aggregate[tool]
            p_ = agg["tp"] / (agg["tp"] + agg["fp"]) if (agg["tp"] + agg["fp"]) else 0.0
            r_ = agg["tp"] / (agg["tp"] + agg["fn"]) if (agg["tp"] + agg["fn"]) else 0.0
            f_ = 2 * p_ * r_ / (p_ + r_) if (p_ + r_) else 0.0
            w.writerow([tool, agg["n_samples"], agg["n_records"],
                        agg["tp"], agg["fp"], agg["fn"],
                        f"{p_:.4f}", f"{r_:.4f}", f"{f_:.4f}"])

    print(f"[score] {len(rows)} per-sample × per-tool rows -> {long_path}")
    print(f"[score] comparative summary -> {comparative_path}")
    # Print summary table
    print()
    print(f"{'tool':<12} {'samples':>7} {'TP':>4} {'FP':>4} {'FN':>4} "
          f"{'precision':>9} {'recall':>7} {'F1':>5}")
    for tool in order:
        agg = per_tool_aggregate[tool]
        p_ = agg["tp"] / (agg["tp"] + agg["fp"]) if (agg["tp"] + agg["fp"]) else 0.0
        r_ = agg["tp"] / (agg["tp"] + agg["fn"]) if (agg["tp"] + agg["fn"]) else 0.0
        f_ = 2 * p_ * r_ / (p_ + r_) if (p_ + r_) else 0.0
        print(f"{tool:<12} {agg['n_samples']:>7} {agg['tp']:>4} {agg['fp']:>4} "
              f"{agg['fn']:>4} {p_:>9.4f} {r_:>7.4f} {f_:>5.4f}")


if __name__ == "__main__":
    main()
