"""quasarsv command-line interface.

Subcommands
-----------
  parse     parse caller VCFs and write per-sample unified TSV
  merge     merge per-caller calls into FusionCall TSV (with tiering)
  annotate  add gene / driver / known-partner annotation
  qc        apply post-merge QC flags (artefact hotspot detection)
  report    render brochure / cohort / validation HTML
  run       end-to-end: VCFs -> merged + annotated + brochure for a sample
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .annotate import annotate_calls
from .merge import MergeConfig, TierThresholds, merge_caller_calls
from .model import (
    FusionCall, BreakpointCall,
    read_fusion_calls_tsv, write_fusion_calls_tsv, write_fusion_calls_json,
)
from .parsers import PARSERS, parse_any
from .qc import apply_default_qc


def _parse_paths(args) -> dict[str, list[str]]:
    """Return {caller: [paths]} from --manta/--gridss/--delly/--svaba/--factera."""
    out: dict[str, list[str]] = {}
    for caller in PARSERS:
        v = getattr(args, caller, None)
        if v:
            out[caller] = list(v)
    return out


def _parse_sample_vcfs(sample: str, paths_by_caller: dict[str, list[str]]) -> list[BreakpointCall]:
    out: list[BreakpointCall] = []
    for caller, paths in paths_by_caller.items():
        for p in paths:
            out.extend(parse_any(p, caller, sample))
    return out


# ---- commands ----

def cmd_parse(args) -> int:
    """Per-caller VCF -> BreakpointCall summary TSV (one row per caller call)."""
    paths_by_caller = _parse_paths(args)
    bps = _parse_sample_vcfs(args.sample, paths_by_caller)
    # Write a flat per-call dump (mostly for debugging / inspection)
    import csv
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "caller", "chrom_a", "pos_a", "strand_a",
            "chrom_b", "pos_b", "strand_b", "sv_type",
            "split_reads", "discordant_pairs", "assembly_contigs",
            "filter_pass", "precise", "raw_qual", "record_id"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for b in bps:
            e = b.evidence
            w.writerow([b.sample, e.caller, b.chrom_a, b.pos_a, b.strand_a,
                        b.chrom_b, b.pos_b, b.strand_b, b.sv_type,
                        e.split_reads, e.discordant_pairs, e.assembly_contigs,
                        e.filter_pass, e.precise, e.raw_qual, b.record_id])
    print(f"[parse] {len(bps)} caller calls -> {out}", file=sys.stderr)
    return 0


def cmd_merge(args) -> int:
    paths_by_caller = _parse_paths(args)
    bps = _parse_sample_vcfs(args.sample, paths_by_caller)
    per_caller: dict[str, list[BreakpointCall]] = {}
    for b in bps:
        per_caller.setdefault(b.evidence.caller, []).append(b)
    cfg = MergeConfig(pos_tolerance=args.tolerance)
    calls = merge_caller_calls(per_caller, cfg=cfg)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_fusion_calls_tsv(calls, str(out))
    print(f"[merge] {len(calls)} candidates -> {out}", file=sys.stderr)
    return 0


def cmd_annotate(args) -> int:
    calls = read_fusion_calls_tsv(args.input)
    annotate_calls(calls)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_fusion_calls_tsv(calls, str(out))
    if args.json:
        write_fusion_calls_json(calls, args.json)
    print(f"[annotate] {len(calls)} annotated -> {out}", file=sys.stderr)
    return 0


def _build_sample_lineage(args) -> tuple[dict[str, str] | None, str]:
    """Return (per-sample lineage map, default lineage) from CLI args.

    Sources, in priority order:
      1. ``--metadata`` XLSX/TSV (per-sample lineage derived from cohort label).
      2. ``--lineage``  fallback default for samples not in the metadata file.
    """
    default = getattr(args, "lineage", None) or "B"
    map_: dict[str, str] | None = None
    meta_path = getattr(args, "metadata", None)
    if meta_path:
        from .metadata import load_cohort_metadata_xlsx, lineage_index_from_metadata
        try:
            items = load_cohort_metadata_xlsx(meta_path)
            map_ = lineage_index_from_metadata(items)
        except Exception as exc:    # noqa: BLE001
            print(f"[lineage] failed to load --metadata {meta_path}: {exc}", file=sys.stderr)
            map_ = None
    return map_, default


def cmd_qc(args) -> int:
    calls = read_fusion_calls_tsv(args.input)
    sample_lineage, lineage_default = _build_sample_lineage(args)
    apply_default_qc(calls, sample_lineage=sample_lineage, lineage_default=lineage_default)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_fusion_calls_tsv(calls, str(out))
    n_flagged = sum(1 for c in calls if c.qc_flags)
    print(f"[qc] {n_flagged}/{len(calls)} flagged -> {out}", file=sys.stderr)
    return 0


def cmd_report(args) -> int:
    from .reports import write_brochure, write_cohort_dashboard, write_validation_report
    calls = read_fusion_calls_tsv(args.input)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    if args.kind in ("brochure", "all"):
        samples = sorted({c.sample for c in calls})
        for s in samples:
            sample_calls = [c for c in calls if c.sample == s]
            p = out_dir / f"brochure_{s}.html"
            write_brochure(s, sample_calls, str(p))
            written.append(str(p))
    if args.kind in ("cohort", "all"):
        p = out_dir / "cohort_dashboard.html"
        write_cohort_dashboard(calls, str(p))
        written.append(str(p))
    if args.kind in ("validation", "all"):
        replicate_pairs = _load_replicate_pairs(args.replicate_pairs)
        p = out_dir / "validation_report.html"
        write_validation_report(calls, str(p), replicate_pairs=replicate_pairs)
        written.append(str(p))
    print(f"[report] wrote {len(written)} reports:", file=sys.stderr)
    for w in written:
        print(f"  {w}", file=sys.stderr)
    return 0


def _load_replicate_pairs(path: str | None) -> list[tuple[str, str]]:
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    pairs: list[tuple[str, str]] = []
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(",")
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    return pairs


def _require_files(*items: tuple[str, str]) -> int:
    """Return non-zero and print a clean message if any (path, label) is missing."""
    for path, label in items:
        if not path or not Path(path).exists():
            print(f"[error] {label} not found: {path}", file=sys.stderr)
            return 2
    return 0


def cmd_scan_cram(args) -> int:
    """Read-level scan of a BAM/CRAM -> merged + annotated TSV + reports."""
    rc = _require_files((args.bam, "BAM/CRAM"), (args.reference, "reference FASTA"))
    if rc:
        return rc
    if not Path(str(args.reference) + ".fai").exists():
        print(f"[error] reference FASTA index missing: {args.reference}.fai "
              f"(run: samtools faidx {args.reference})", file=sys.stderr)
        return 2
    from .scanners import ScannerConfig, SAScannerConfig, scan_cram, scan_artefacts_sa
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = ScannerConfig(
        min_mapq=args.min_mapq,
        min_split_reads=args.min_split_reads,
        min_discordant_pairs=args.min_discordant_pairs,
        pad_locus_bp=args.pad_locus_bp,
        use_mapq_weighting=not args.no_mapq_weighting,
        use_adaptive_insert=not args.no_adaptive_insert,
        full_weight_mapq=args.full_weight_mapq,
    )
    lib_path = (str(out_dir / "library_stats.json")
                if not args.no_adaptive_insert else None)
    bps = scan_cram(args.bam, args.reference, args.sample, cfg=cfg,
                    library_stats_path=lib_path)
    # SA-aware artefact scan — genuine chimeras that overlap an artefact window.
    # Poly-G/adapter alignments inside the attractor are filtered out; they name
    # the junk read's origin locus, not a translocation partner.
    sa_bps = scan_artefacts_sa(args.bam, args.reference, args.sample,
                               cfg=SAScannerConfig(min_split_reads=max(5, args.min_split_reads)))
    per_caller = {"forge_scan": bps, "forge_scan_sa": sa_bps}
    chrom_sa_bps: list = []
    calls = merge_caller_calls(per_caller, cfg=MergeConfig(pos_tolerance=250))
    annotate_calls(calls)
    sample_lineage, lineage_default = _build_sample_lineage(args)
    apply_default_qc(calls, sample_lineage=sample_lineage, lineage_default=lineage_default)
    tsv = out_dir / f"{args.sample}.fusions.tsv"
    write_fusion_calls_tsv(calls, str(tsv))
    write_fusion_calls_json(calls, str(out_dir / f"{args.sample}.fusions.json"))

    summary = {
        "sample": args.sample,
        "n_breakpoint_calls": len(bps) + len(sa_bps) + len(chrom_sa_bps),
        "n_forge_scan": len(bps),
        "n_forge_scan_sa": len(sa_bps),
        "n_forge_scan_chrom_sa": len(chrom_sa_bps),
        "n_candidates": len(calls),
        "tier": {t: sum(1 for c in calls if c.tier == t) for t in ("T1", "T2", "T3")},
        "known_partner_count": sum(1 for c in calls if c.known_partner),
        "driver_T1_T2": sum(1 for c in calls if c.driver_locus and c.tier in ("T1", "T2")),
    }
    if not args.skip_reports:
        from .reports import write_brochure, write_cohort_dashboard, write_validation_report
        brochure = out_dir / f"brochure_{args.sample}.html"
        write_brochure(args.sample, calls, str(brochure))
        write_cohort_dashboard(calls, str(out_dir / "cohort_dashboard.html"))
        write_validation_report(calls, str(out_dir / "validation_report.html"))
        summary["brochure"] = str(brochure)
    print(json.dumps(summary, indent=2))
    return 0


def cmd_emit_vcf(args) -> int:
    """Emit quasarsv FusionCalls as a VCF 4.3 file (gz-aware via extension)."""
    from .vcf_emit import write_vcf
    calls = read_fusion_calls_tsv(args.input)
    if not calls:
        print(f"[emit-vcf] no calls in {args.input}", file=sys.stderr)
        return 1
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    n = write_vcf(calls, str(out), sample=args.sample, emit_mates=not args.no_mates)
    print(f"[emit-vcf] {n} VCF records -> {out}", file=sys.stderr)
    return 0


def cmd_benchmark(args) -> int:
    """Score one or more FusionCall TSVs against a truth-set TSV."""
    from .benchmark import (
        builtin_truth_path, load_truth_set, score_calls_against_truth,
        write_benchmark_tsv,
    )
    truth_path = args.truth_set or str(builtin_truth_path())
    truth = load_truth_set(truth_path)
    if not truth:
        print(f"[benchmark] truth set empty / not found: {truth_path}", file=sys.stderr)
        return 1

    calls: list[FusionCall] = []
    for p in args.inputs:
        calls.extend(read_fusion_calls_tsv(p))
    if not calls:
        print(f"[benchmark] no calls loaded from {args.inputs}", file=sys.stderr)
        return 1

    junction_support = None
    if getattr(args, "junction_truth", None):
        from .benchmark import load_junction_support
        junction_support = load_junction_support(args.junction_truth)

    bench = score_calls_against_truth(
        calls, truth,
        tool_name=args.tool_name,
        ambiguous_alts_count_as_fp=args.ambiguous_alts_count_as_fp,
        relax_canonical_ig_partner=args.relax_canonical_ig_partner,
        match_driver_only=getattr(args, "match_driver_only", False),
        junction_support=junction_support,
        restrict_to_scored_samples=not args.include_unscanned,
    )
    out_path = Path(args.output) if args.output else None
    if out_path:
        write_benchmark_tsv(bench, out_path)

    samples_with_truth = len(bench.per_sample)
    samples_scored = sum(1 for s in bench.per_sample.values() if s.expected_pairs)
    samples_neg = sum(1 for s in bench.per_sample.values() if s.is_negative_control)
    print(json.dumps({
        "tool": bench.tool_name,
        "samples_in_truth": samples_with_truth,
        "samples_with_positive_truth": samples_scored,
        "samples_negative_control": samples_neg,
        "tp": bench.tp, "fp": bench.fp, "fn": bench.fn,
        "precision": round(bench.precision, 4),
        "recall": round(bench.recall, 4),
        "f1": round(bench.f1, 4),
        "output_tsv": str(out_path) if out_path else None,
    }, indent=2))
    return 0


def cmd_run(args) -> int:
    """End-to-end: parse → merge → annotate → qc → reports."""
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths_by_caller = _parse_paths(args)
    bps = _parse_sample_vcfs(args.sample, paths_by_caller)
    per_caller: dict[str, list[BreakpointCall]] = {}
    for b in bps:
        per_caller.setdefault(b.evidence.caller, []).append(b)
    calls = merge_caller_calls(per_caller, cfg=MergeConfig(pos_tolerance=args.tolerance))
    annotate_calls(calls)
    sample_lineage, lineage_default = _build_sample_lineage(args)
    apply_default_qc(calls, sample_lineage=sample_lineage, lineage_default=lineage_default)

    tsv = out_dir / f"{args.sample}.fusions.tsv"
    write_fusion_calls_tsv(calls, str(tsv))
    js = out_dir / f"{args.sample}.fusions.json"
    write_fusion_calls_json(calls, str(js))

    from .reports import write_brochure, write_cohort_dashboard, write_validation_report
    brochure = out_dir / f"brochure_{args.sample}.html"
    write_brochure(args.sample, calls, str(brochure))
    cohort = out_dir / f"cohort_dashboard.html"
    write_cohort_dashboard(calls, str(cohort))
    val = out_dir / f"validation_report.html"
    write_validation_report(calls, str(val))

    summary = {
        "sample": args.sample,
        "n_caller_records": len(bps),
        "n_candidates": len(calls),
        "tier": {t: sum(1 for c in calls if c.tier == t) for t in ("T1", "T2", "T3")},
        "known_partner_count": sum(1 for c in calls if c.known_partner),
        "outputs": {"tsv": str(tsv), "json": str(js),
                     "brochure": str(brochure), "cohort": str(cohort), "validation": str(val)},
    }
    print(json.dumps(summary, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="quasar",
        description="Quasar — lymphoma-specific structural-variant / fusion caller "
                    "(BAM/CRAM in; VCF + tiered TSV + HTML report out)")
    sub = parser.add_subparsers(dest="cmd", required=True)

    def add_caller_args(p: argparse.ArgumentParser):
        for c in PARSERS:
            p.add_argument(f"--{c}", nargs="*", help=f"{c} VCF / output paths")

    p_parse = sub.add_parser("parse", help="dump per-caller calls to TSV")
    add_caller_args(p_parse)
    p_parse.add_argument("--sample", required=True)
    p_parse.add_argument("--output", required=True)
    p_parse.set_defaults(func=cmd_parse)

    p_merge = sub.add_parser("merge", help="merge caller VCFs into FusionCall TSV")
    add_caller_args(p_merge)
    p_merge.add_argument("--sample", required=True)
    p_merge.add_argument("--output", required=True)
    p_merge.add_argument("--tolerance", type=int, default=250,
                         help="bp window for clustering breakpoints across callers")
    p_merge.set_defaults(func=cmd_merge)

    p_ann = sub.add_parser("annotate", help="annotate a FusionCall TSV with genes / known partners")
    p_ann.add_argument("--input", required=True)
    p_ann.add_argument("--output", required=True)
    p_ann.add_argument("--json", help="optional JSON sibling output")
    p_ann.set_defaults(func=cmd_annotate)

    def _add_lineage_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("--lineage", choices=["B", "T", "any"], default="B",
                       help="default lymphoma lineage prior for the artefact-rescue "
                            "IG partner selection (B=IGH/IGK/IGL only, T=TRA/TRB/TRG/TRD only; "
                            "default B)")
        p.add_argument("--metadata",
                       help="cohort metadata XLSX with a 'Cohort' column to auto-infer "
                            "per-sample lineage (overrides --lineage per sample)")

    p_qc = sub.add_parser("qc", help="apply post-merge QC flags")
    p_qc.add_argument("--input", required=True)
    p_qc.add_argument("--output", required=True)
    _add_lineage_args(p_qc)
    p_qc.set_defaults(func=cmd_qc)

    p_rep = sub.add_parser("report", help="render HTML reports from a FusionCall TSV")
    p_rep.add_argument("--input", required=True)
    p_rep.add_argument("--output-dir", required=True)
    p_rep.add_argument("--kind", choices=["brochure", "cohort", "validation", "all"], default="all")
    p_rep.add_argument("--replicate-pairs",
                       help="TSV/CSV of replicate sample pairs (one pair per line)")
    p_rep.set_defaults(func=cmd_report)

    p_vcf = sub.add_parser("emit-vcf",
                            help="convert a FusionCall TSV to a VCF 4.3 file")
    p_vcf.add_argument("--input", required=True, help="*.fusions.tsv path")
    p_vcf.add_argument("--output", required=True, help="*.vcf or *.vcf.gz path")
    p_vcf.add_argument("--sample", help="sample id for FORMAT/sample column")
    p_vcf.add_argument("--no-mates", action="store_true",
                        help="emit one record per BND instead of both mate ends")
    p_vcf.set_defaults(func=cmd_emit_vcf)

    p_bench = sub.add_parser("benchmark",
                              help="score a FusionCall TSV against a truth-set TSV")
    p_bench.add_argument("inputs", nargs="+",
                          help="one or more *.fusions.tsv to score")
    p_bench.add_argument("--truth-set",
                          help="truth-set TSV path (defaults to packaged data/cohort_truth.tsv)")
    p_bench.add_argument("--tool-name", default="quasarsv",
                          help="label written into the benchmark output (e.g. manta, gridss2)")
    p_bench.add_argument("--output", help="output benchmark TSV (per-sample + aggregate)")
    p_bench.add_argument("--no-ambiguous-alts-count-as-fp",
                          dest="ambiguous_alts_count_as_fp",
                          action="store_false", default=True,
                          help="exempt ig_partner_ambiguous T1 calls from false positives; "
                               "default OFF — hedged alternates count against precision")
    p_bench.add_argument("--relax-canonical-ig-partner", action="store_true",
                          help="LENIENT side-metric, never the headline: BCL6-IGL matches "
                               "truth BCL6-IGH when both IGs are canonical. Measures "
                               "driver identification, not partner resolution")
    p_bench.add_argument("--match-driver-only", action="store_true",
                          help="LENIENT side-metric: a driver-only call matches a driver-IG "
                               "truth. Requires --relax-canonical-ig-partner")
    p_bench.add_argument("--junction-truth",
                          help="TSV of independently-established read-level junctions "
                               "(sample, gene_a, gene_b, ...) used to split true positives "
                               "into detected vs lookup-only")
    p_bench.add_argument("--include-unscanned", action="store_true",
                          help="include truth-set samples with no calls (counted as FN); "
                               "default off — unscanned ≠ missed")
    p_bench.set_defaults(func=cmd_benchmark)

    p_run = sub.add_parser("run", help="end-to-end VCFs -> reports for a single sample")
    add_caller_args(p_run)
    p_run.add_argument("--sample", required=True)
    p_run.add_argument("--output-dir", required=True)
    p_run.add_argument("--tolerance", type=int, default=250)
    _add_lineage_args(p_run)
    p_run.set_defaults(func=cmd_run)

    p_scan = sub.add_parser("call",
                            aliases=["scan-cram"],
                            help="call structural variants from a BAM/CRAM")
    p_scan.add_argument("--sample", required=True)
    p_scan.add_argument("--bam", required=True, help="BAM or CRAM path")
    p_scan.add_argument("--reference", required=True, help="FASTA matching the CRAM @SQ")
    p_scan.add_argument("--output-dir", required=True)
    p_scan.add_argument("--min-mapq", type=int, default=20)
    p_scan.add_argument("--min-split-reads", type=int, default=2)
    p_scan.add_argument("--min-discordant-pairs", type=int, default=4)
    p_scan.add_argument("--pad-locus-bp", type=int, default=5_000)
    p_scan.add_argument("--skip-reports", action="store_true",
                        help="only emit the fusion TSV; skip HTML render")
    p_scan.add_argument("--chrom-sa-inference", action="store_true",
                        help=argparse.SUPPRESS)   # removed: fabricated breakpoint coordinates
    p_scan.add_argument("--no-mapq-weighting", action="store_true",
                        help="disable GRIDSS-style MAPQ-weighted contribution; revert to "
                             "hard-cutoff counting (legacy mode)")
    p_scan.add_argument("--no-adaptive-insert", action="store_true",
                        help="disable Delly-style adaptive insert threshold inferred from "
                             "library stats; revert to fixed discordant_min_distance")
    p_scan.add_argument("--full-weight-mapq", type=int, default=30,
                        help="MAPQ at and above which a read contributes weight 1.0 "
                             "(below tapers linearly to 0); default 30")
    _add_lineage_args(p_scan)
    p_scan.set_defaults(func=cmd_scan_cram)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
