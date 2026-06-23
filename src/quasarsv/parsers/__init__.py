"""Caller-specific VCF parsers — each emits BreakpointCall in the unified model."""
from .base import open_vcf, parse_bnd_alt, normalise_order
from .manta import parse_manta
from .gridss import parse_gridss
from .delly import parse_delly
from .svaba import parse_svaba
from .factera import parse_factera
from .tiddit import parse_tiddit

PARSERS = {
    "manta": parse_manta,
    "gridss": parse_gridss,
    "delly": parse_delly,
    "svaba": parse_svaba,
    "factera": parse_factera,
    "tiddit": parse_tiddit,
}


def parse_any(path: str, caller: str, sample: str):
    if caller not in PARSERS:
        raise ValueError(f"Unknown caller {caller!r}; expected one of {list(PARSERS)}")
    return PARSERS[caller](path, sample)
