"""Shared report styling and HTML helpers — dark navy header, card sections, tables-first."""
from __future__ import annotations

import html
from datetime import datetime, timezone


BROCHURE_CSS = """
:root{
  --bg:#fafafa; --fg:#1a1a1a; --muted:#666; --accent:#d62728;
  --card:#fff; --navy:#16213e; --rule:#e8e8e8;
}
*{box-sizing:border-box}
body{font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;
     background:var(--bg);color:var(--fg);margin:0;line-height:1.5}
header{background:linear-gradient(135deg,#1a1a2e,#16213e);color:white;padding:36px 32px}
header h1{margin:0;font-size:30px;letter-spacing:-0.01em}
header p{margin:6px 0 0;opacity:0.82;font-size:13px}
nav{background:#fff;padding:10px 32px;border-bottom:1px solid var(--rule);
    position:sticky;top:0;z-index:10;font-size:13px}
nav a{color:var(--navy);text-decoration:none;margin-right:18px;font-weight:500}
nav a:hover{color:var(--accent)}
main{max-width:1220px;margin:0 auto;padding:24px}
.section{background:var(--card);border-radius:8px;
         box-shadow:0 1px 3px rgba(0,0,0,0.06);padding:24px;margin-bottom:20px}
.section h2{margin:0 0 10px;font-size:20px;color:var(--navy);
            border-bottom:2px solid var(--accent);padding-bottom:6px;display:inline-block}
.section h3{margin:18px 0 8px;font-size:15px;color:var(--navy)}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
      gap:14px;margin:14px 0}
.kpi{background:#fff;border:1px solid var(--rule);border-radius:6px;
     padding:14px;text-align:center}
.kpi .v{font-size:28px;font-weight:600;color:var(--navy);font-variant-numeric:tabular-nums}
.kpi .l{font-size:11px;color:var(--muted);margin-top:4px;
        text-transform:uppercase;letter-spacing:0.5px}
table.dt{width:100%;border-collapse:collapse;margin-top:10px;font-size:13px;
         font-variant-numeric:tabular-nums}
table.dt th{background:var(--navy);color:white;padding:8px 10px;text-align:left;
            font-weight:600}
table.dt td{padding:6px 10px;border-bottom:1px solid var(--rule);vertical-align:top}
table.dt tbody tr:hover{background:#f0f4f8}
.tier{display:inline-block;padding:2px 8px;border-radius:10px;font-weight:600;
      font-size:11px;color:white;letter-spacing:0.4px}
.tier.T1{background:#d62728}.tier.T2{background:#ff7f0e}.tier.T3{background:#9ecae1;color:#333}
.badge{display:inline-block;padding:2px 6px;font-size:11px;border-radius:4px;
       background:#eef2f6;color:var(--navy);margin-right:4px;font-weight:600;letter-spacing:0.3px}
.badge.warn{background:#fff3cd;color:#7a5a00}
.badge.partner{background:#dbeafe;color:#0c4a6e}
.badge.ig{background:#0c4a6e;color:white}
.badge.tr{background:#5b21b6;color:white}
.badge.drv{background:#b91c1c;color:white}
.translocation-card{border:1px solid var(--rule);border-radius:6px;
                    padding:14px;margin-top:14px;background:#fcfcfd}
.translocation-card h3{margin:0 0 6px;font-size:15px;color:var(--navy);display:flex;
                       align-items:center;gap:8px}
.translocation-card .cyto{font-family:'SFMono-Regular',Menlo,Consolas,monospace;
                          font-size:12px;color:var(--muted)}
.translocation-card .disease{font-size:12px;color:var(--muted);margin-bottom:4px}
.empty-row{color:var(--muted);font-size:12px;padding:4px 0}
.plot-wrap{margin-top:12px;border:1px solid var(--rule);border-radius:6px;
           padding:6px;background:#fff}
footer{text-align:center;color:var(--muted);padding:22px;font-size:12px}
.small{font-size:12px;color:var(--muted)}
.mono{font-family:'SFMono-Regular',Menlo,Consolas,monospace;font-size:12px}
"""

PLOTLY_JS = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'


def now_str() -> str:
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M %Z")


def fig_to_div(fig, div_id: str | None = None) -> str:
    """Embed a plotly figure as a self-contained <div>."""
    import plotly.io as pio
    return pio.to_html(fig, include_plotlyjs=False, full_html=False,
                       div_id=div_id, default_width="100%", default_height=None,
                       config={"displayModeBar": False, "responsive": True})


def render_table(rows: list[dict], columns: list[tuple[str, str]],
                 row_class: callable | None = None) -> str:
    """columns is list of (column_label, dict_key)."""
    if not rows:
        return "<p class='small'>No rows.</p>"
    th = "".join(f"<th>{html.escape(lbl)}</th>" for lbl, _ in columns)
    body_rows = []
    for r in rows:
        cls = row_class(r) if row_class else ""
        tds = "".join(f"<td>{_fmt(r.get(k))}</td>" for _, k in columns)
        body_rows.append(f"<tr class='{cls}'>{tds}</tr>")
    return f"<table class='dt'><thead><tr>{th}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


class Safe(str):
    """Marker for strings that already contain trusted HTML — bypasses escaping."""
    __slots__ = ()


def _fmt(v) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, Safe):
        return str(v)
    if isinstance(v, bool):
        return "✓" if v else ""
    if isinstance(v, float):
        if abs(v) < 1e-4 and v != 0:
            return f"{v:.2e}"
        return f"{v:.3f}".rstrip("0").rstrip(".")
    if isinstance(v, (list, tuple)):
        return ", ".join(html.escape(str(x)) for x in v)
    return html.escape(str(v))


def tier_badge(tier: str) -> Safe:
    return Safe(f"<span class='tier {tier}'>{tier}</span>")


_IG_SET = {"IGH", "IGK", "IGL", "IGH_Emu", "IGH_3RR"}
_TR_SET = {"TRA", "TRB", "TRG", "TRD"}


def locus_badge(gene: str) -> Safe:
    """Coloured chip for a gene, distinguishing IG / TR / driver / other."""
    if not gene:
        return Safe("<span class='badge'>intergenic</span>")
    if gene in _IG_SET:
        return Safe(f"<span class='badge ig'>IG · {gene}</span>")
    if gene in _TR_SET:
        return Safe(f"<span class='badge tr'>TR · {gene}</span>")
    return Safe(f"<span class='badge drv'>{gene}</span>")


def gene_pair_with_badges(gene_a: str, gene_b: str) -> Safe:
    """Both endpoints with IG / TR / driver badge plus ↔ between them."""
    return Safe(f"{locus_badge(gene_a)} ↔ {locus_badge(gene_b)}")


def kpi(label: str, value) -> str:
    return f"<div class='kpi'><div class='v'>{html.escape(str(value))}</div><div class='l'>{html.escape(label)}</div></div>"


def kpi_grid(items: list[tuple[str, str]]) -> str:
    return "<div class='kpis'>" + "".join(kpi(l, v) for l, v in items) + "</div>"


def html_shell(title: str, nav_items: list[tuple[str, str]], body: str) -> str:
    nav = "".join(f"<a href='#{aid}'>{html.escape(lbl)}</a>" for lbl, aid in nav_items)
    return f"""<!doctype html><html><head><meta charset='utf-8'>
<title>{html.escape(title)}</title>
{PLOTLY_JS}
<style>{BROCHURE_CSS}</style>
</head><body>
<header><h1>{html.escape(title)}</h1><p>generated {now_str()}</p></header>
<nav>{nav}</nav>
<main>{body}</main>
<footer>quasarsv v0.1.0 — built for lymphoma fusion / rearrangement detection from NGS data</footer>
</body></html>"""
