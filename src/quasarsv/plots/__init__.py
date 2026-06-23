"""Plot suite. Each module returns a plotly Figure or HTML div string."""
from .circos import circos_figure
from .validation import per_caller_concordance, per_caller_evidence_breakdown, replicate_concordance
from .qc_plots import qc_summary_figure, evidence_distribution
from .locus import locus_figure, locus_summary_table

__all__ = [
    "circos_figure",
    "per_caller_concordance",
    "per_caller_evidence_breakdown",
    "replicate_concordance",
    "qc_summary_figure",
    "evidence_distribution",
    "locus_figure",
    "locus_summary_table",
]
