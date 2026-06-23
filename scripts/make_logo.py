#!/usr/bin/env python3
"""Build the Quasar visual identity.

Outputs under docs/branding/:
  quasar_icon.svg / .png            square icon, 1024x1024 PNG
  quasar_logo.svg / .png            wordmark + tagline, light theme
  quasar_logo_dark.png              wordmark + tagline, dark theme
  quasar_favicon.ico                multi-resolution (16, 32, 48, 64)

The SVGs use real radialGradient / linearGradient and are the master assets.
PNGs are rendered via Pillow with manual radial-blend compositing.
"""
from __future__ import annotations

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "branding"
OUT.mkdir(parents=True, exist_ok=True)

# Palette
NAVY      = (11, 30, 58)        # #0b1e3a — light-theme text + accents
DEEP      = (5, 10, 24)         # #050a18 — dark-theme background
GOLD_HOT  = (255, 244, 200)     # core hottest
GOLD      = (251, 191, 36)      # core mid
GOLD_DIM  = (217, 119, 6)       # core outer
DISK_OUTER = (124, 58, 237)     # #7c3aed violet
DISK_INNER = (236, 72, 153)     # #ec4899 magenta
JET_HOT   = (255, 255, 255)     # jet centre
JET_COOL  = (34, 211, 238)      # #22d3ee cyan
SUB_LIGHT = (71, 85, 105)       # #475569 subtitle on light
SUB_DARK  = (148, 163, 184)     # #94a3b8 subtitle on dark


# ---------------------------------------------------------------------------
# SVG master (hand-written, gradient-rich)
# ---------------------------------------------------------------------------

ICON_SVG = '''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">
  <defs>
    <!-- bright core -->
    <radialGradient id="core" cx="50%" cy="50%" r="50%">
      <stop offset="0%"  stop-color="#fff7d6" stop-opacity="1"/>
      <stop offset="35%" stop-color="#fde68a" stop-opacity="0.95"/>
      <stop offset="65%" stop-color="#fbbf24" stop-opacity="0.55"/>
      <stop offset="100%" stop-color="#fbbf24" stop-opacity="0"/>
    </radialGradient>
    <!-- accretion disk: violet -> magenta sweep -->
    <linearGradient id="disk" x1="0%" y1="50%" x2="100%" y2="50%">
      <stop offset="0%"   stop-color="#7c3aed"/>
      <stop offset="50%"  stop-color="#ec4899"/>
      <stop offset="100%" stop-color="#7c3aed"/>
    </linearGradient>
    <!-- relativistic jet -->
    <linearGradient id="jet" x1="50%" y1="0%" x2="50%" y2="100%">
      <stop offset="0%"   stop-color="#22d3ee" stop-opacity="0"/>
      <stop offset="20%"  stop-color="#22d3ee" stop-opacity="0.65"/>
      <stop offset="50%"  stop-color="#ffffff" stop-opacity="0.95"/>
      <stop offset="80%"  stop-color="#22d3ee" stop-opacity="0.65"/>
      <stop offset="100%" stop-color="#22d3ee" stop-opacity="0"/>
    </linearGradient>
    <!-- soft outer halo -->
    <radialGradient id="halo" cx="50%" cy="50%" r="50%">
      <stop offset="0%"   stop-color="#fbbf24" stop-opacity="0.35"/>
      <stop offset="60%"  stop-color="#fbbf24" stop-opacity="0.05"/>
      <stop offset="100%" stop-color="#fbbf24" stop-opacity="0"/>
    </radialGradient>
  </defs>
  <!-- transparent background -->
  <rect width="1024" height="1024" fill="none"/>
  <!-- soft golden halo behind everything -->
  <circle cx="512" cy="512" r="280" fill="url(#halo)"/>
  <!-- jet: tapered bowtie via polygon, filled with vertical linear gradient -->
  <polygon points="430,40 594,40 542,512 482,512" fill="url(#jet)"/>
  <polygon points="430,984 594,984 542,512 482,512" fill="url(#jet)"/>
  <!-- accretion disk: ellipse stroke + glow ring -->
  <ellipse cx="512" cy="512" rx="380" ry="64"
           fill="none" stroke="url(#disk)" stroke-width="14" opacity="0.95"/>
  <ellipse cx="512" cy="512" rx="380" ry="64"
           fill="none" stroke="url(#disk)" stroke-width="40" opacity="0.18"/>
  <ellipse cx="512" cy="512" rx="260" ry="38"
           fill="none" stroke="url(#disk)" stroke-width="6"  opacity="0.7"/>
  <!-- central engine: glowing core -->
  <circle cx="512" cy="512" r="140" fill="url(#core)"/>
  <circle cx="512" cy="512" r="22"  fill="#ffffff"/>
</svg>
'''


def write_svg(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    print(f"[svg]   {path}")


# ---------------------------------------------------------------------------
# Pillow renderers
# ---------------------------------------------------------------------------

def _radial_gradient(size: int, stops: list[tuple[float, tuple[int, int, int, int]]]) -> Image.Image:
    """Render an RGBA circle of side ``size`` with a radial gradient.

    stops is a list of (radius_fraction_0_to_1, (R, G, B, A))."""
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    cx = cy = size / 2.0
    max_r = size / 2.0
    px = img.load()
    for y in range(size):
        for x in range(size):
            dx, dy = x - cx, y - cy
            r = math.sqrt(dx * dx + dy * dy) / max_r
            if r > 1.0:
                continue
            # find bracketing stops + lerp
            prev = stops[0]
            for stop in stops[1:]:
                if r <= stop[0]:
                    span = stop[0] - prev[0]
                    t = 0.0 if span == 0 else (r - prev[0]) / span
                    c0 = prev[1]; c1 = stop[1]
                    px[x, y] = (
                        int(c0[0] * (1 - t) + c1[0] * t),
                        int(c0[1] * (1 - t) + c1[1] * t),
                        int(c0[2] * (1 - t) + c1[2] * t),
                        int(c0[3] * (1 - t) + c1[3] * t),
                    )
                    break
                prev = stop
    return img


def _draw_jet(canvas: Image.Image, cx: int, cy: int, length: int,
              half_width_base: int, half_width_tip: int) -> None:
    """Two tapered jets, top and bottom, with cyan→white→cyan vertical gradient."""
    for direction in (1, -1):
        # Build a thin vertical gradient strip then warp via polygon mask
        strip_w = max(half_width_base * 2, 16)
        strip_h = length
        strip = Image.new("RGBA", (strip_w, strip_h), (0, 0, 0, 0))
        sp = strip.load()
        for y in range(strip_h):
            t = y / max(strip_h - 1, 1)
            # cyan -> white -> cyan along axis; alpha fades toward ends
            if t < 0.5:
                a = t / 0.5
                r = int(JET_COOL[0] * (1 - a) + JET_HOT[0] * a)
                g = int(JET_COOL[1] * (1 - a) + JET_HOT[1] * a)
                b = int(JET_COOL[2] * (1 - a) + JET_HOT[2] * a)
            else:
                a = (t - 0.5) / 0.5
                r = int(JET_HOT[0] * (1 - a) + JET_COOL[0] * a)
                g = int(JET_HOT[1] * (1 - a) + JET_COOL[1] * a)
                b = int(JET_HOT[2] * (1 - a) + JET_COOL[2] * a)
            # alpha — fades at the far end (away from core)
            alpha_axis = 1.0 - abs(t - 0.5) * 1.4  # peak in middle
            alpha_axis = max(0, min(1, alpha_axis))
            for x in range(strip_w):
                # taper: width shrinks from full at core (y=length-1 if direction=+1
                # i.e. attached to core) toward tip
                if direction == 1:
                    progress = y / max(strip_h - 1, 1)  # 0 at top tip, 1 at core base
                else:
                    progress = 1 - y / max(strip_h - 1, 1)
                half_w = half_width_tip + (half_width_base - half_width_tip) * progress
                if half_w < 0.5:
                    continue
                center = strip_w / 2
                dx = abs(x - center)
                if dx > half_w:
                    continue
                # soft edge
                edge = max(0, 1.0 - (dx / half_w) ** 1.8)
                alpha = int(255 * alpha_axis * edge)
                if alpha > 0:
                    sp[x, y] = (r, g, b, alpha)
        # blit jet onto canvas at correct position
        if direction == 1:
            # top jet: tip at (cx, cy-length), base at (cx, cy)
            canvas.alpha_composite(strip, (cx - strip_w // 2, cy - length))
        else:
            # bottom jet: base at (cx, cy), tip at (cx, cy+length)
            canvas.alpha_composite(strip, (cx - strip_w // 2, cy))


def _draw_disk(draw: ImageDraw.ImageDraw, cx: int, cy: int,
               rx: int, ry: int, stroke_width: int,
               color: tuple[int, int, int, int]) -> None:
    """Stroked ellipse — Pillow has no native stroke-only ellipse, so we draw
    a filled outer ellipse and erase the interior."""
    if rx <= 1 or ry <= 1:
        return
    pad = stroke_width // 2
    inner_rx = rx - pad
    inner_ry = ry - pad
    bbox_outer = [cx - rx - pad, cy - ry - pad, cx + rx + pad, cy + ry + pad]
    draw.ellipse(bbox_outer, fill=color)
    if inner_rx > 0 and inner_ry > 0:
        bbox_inner = [cx - inner_rx, cy - inner_ry, cx + inner_rx, cy + inner_ry]
        draw.ellipse(bbox_inner, fill=(0, 0, 0, 0))


def render_icon_png(size: int = 1024, on_dark: bool = False,
                    compact: bool = False) -> Image.Image:
    """Render the Quasar icon.

    compact=True tightens jets + disk so the icon reads as a single unit
    when placed next to a wordmark.
    """
    if compact:
        s_halo, s_jet, s_jet_base, s_jet_tip = 0.42, 0.50, 0.075, 0.020
        s_disk_x, s_disk_y, s_disk_stroke = 0.46, 0.085, 0.022
        s_core = 0.42
    else:
        s_halo, s_jet, s_jet_base, s_jet_tip = 0.55, 0.46, 0.052, 0.015
        s_disk_x, s_disk_y, s_disk_stroke = 0.37, 0.063, 0.014
        s_core = 0.30

    bg = DEEP + (255,) if on_dark else (255, 255, 255, 0)
    img = Image.new("RGBA", (size, size), bg)
    cx = cy = size // 2

    halo = _radial_gradient(
        int(size * s_halo),
        [(0.0, (251, 191, 36, 90)),
         (0.6, (251, 191, 36, 12)),
         (1.0, (251, 191, 36, 0))],
    )
    img.alpha_composite(halo, (cx - halo.width // 2, cy - halo.height // 2))

    jet_length = int(size * s_jet)
    _draw_jet(img, cx, cy, jet_length,
              half_width_base=max(8, int(size * s_jet_base)),
              half_width_tip=max(3, int(size * s_jet_tip)))

    draw = ImageDraw.Draw(img)
    _draw_disk(draw, cx, cy, int(size * s_disk_x), int(size * s_disk_y),
               stroke_width=int(size * 0.038), color=DISK_OUTER + (45,))
    _draw_disk(draw, cx, cy, int(size * s_disk_x), int(size * s_disk_y),
               stroke_width=int(size * s_disk_stroke), color=DISK_OUTER + (240,))
    _draw_disk(draw, cx, cy, int(size * s_disk_x * 0.69),
               int(size * s_disk_y * 0.58),
               stroke_width=max(2, int(size * 0.008)), color=DISK_INNER + (180,))

    core_size = int(size * s_core)
    core = _radial_gradient(
        core_size,
        [(0.0,  GOLD_HOT + (255,)),
         (0.18, GOLD_HOT + (240,)),
         (0.40, GOLD + (200,)),
         (0.70, GOLD_DIM + (90,)),
         (1.0,  GOLD_DIM + (0,))],
    )
    img.alpha_composite(core, (cx - core.width // 2, cy - core.height // 2))

    pin_r = max(2, int(size * 0.022))
    draw.ellipse([cx - pin_r, cy - pin_r, cx + pin_r, cy + pin_r],
                 fill=(255, 255, 255, 255))

    return img


def _try_font(family: str, size: int) -> ImageFont.FreeTypeFont:
    """Resolve a font by family with cross-platform fallbacks."""
    candidates = [
        family,
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/mnt/c/Windows/Fonts/segoeuib.ttf",
        "/mnt/c/Windows/Fonts/arialbd.ttf",
        "/mnt/c/Windows/Fonts/segoeui.ttf",
        "/mnt/c/Windows/Fonts/arial.ttf",
    ]
    for c in candidates:
        try:
            return ImageFont.truetype(c, size)
        except (OSError, IOError):
            continue
    return ImageFont.load_default()


def render_full_logo(height: int = 360, on_dark: bool = False) -> Image.Image:
    """Horizontal logo: icon on the left, wordmark + tagline on the right.

    Canvas width is sized to fit the actual text content (no clipping).
    Icon height ≈ full logo height; wordmark cap-height ≈ icon height for
    visual balance.
    """
    bg = DEEP if on_dark else (255, 255, 255)
    fg = (255, 255, 255) if on_dark else NAVY
    sub = SUB_DARK if on_dark else SUB_LIGHT

    # icon at full height (compact variant — tighter jets + larger disk)
    icon_size = int(height * 1.0)
    pad_left = int(height * 0.08)
    gap_icon_text = int(height * 0.06)

    # type sizes — wordmark sized so cap height ≈ icon's visual core
    wm_font_size = int(height * 0.42)
    sub_font_size = int(height * 0.095)
    wm_font = _try_font("segoeuib", wm_font_size)
    sub_font = _try_font("segoeui", sub_font_size)

    wm_text = "Quasar"
    sub_text = "lymphoma structural-variant and fusion detection"

    # Measure text on a throwaway canvas to size the real one
    probe = Image.new("RGBA", (10, 10))
    pd = ImageDraw.Draw(probe)
    wm_bbox = pd.textbbox((0, 0), wm_text, font=wm_font)
    sub_bbox = pd.textbbox((0, 0), sub_text, font=sub_font)
    wm_w = wm_bbox[2] - wm_bbox[0]
    sub_w = sub_bbox[2] - sub_bbox[0]
    text_w = max(wm_w, sub_w)

    pad_right = int(height * 0.18)
    width = pad_left + icon_size + gap_icon_text + text_w + pad_right

    img = Image.new("RGBA", (width, height), bg + (255,))

    # icon — use compact variant for better visual balance vs wordmark
    icon = render_icon_png(size=icon_size, on_dark=on_dark, compact=True)
    icon_y = (height - icon_size) // 2
    img.alpha_composite(icon, (pad_left, icon_y))

    draw = ImageDraw.Draw(img)
    text_x = pad_left + icon_size + gap_icon_text

    wm_h = wm_bbox[3] - wm_bbox[1]
    sub_h = sub_bbox[3] - sub_bbox[1]
    gap_text = int(height * 0.04)
    block_h = wm_h + gap_text + sub_h
    wm_y = (height - block_h) // 2 - wm_bbox[1]
    sub_y = wm_y + wm_h + gap_text - sub_bbox[1]

    draw.text((text_x, wm_y), wm_text, font=wm_font, fill=fg + (255,))
    draw.text((text_x, sub_y), sub_text, font=sub_font, fill=sub + (255,))

    return img


def main() -> None:
    # Icon SVG + PNG (transparent bg)
    write_svg(OUT / "quasar_icon.svg", ICON_SVG)

    icon = render_icon_png(size=1024, on_dark=False)
    icon.save(OUT / "quasar_icon.png", "PNG")
    print(f"[png]   {OUT / 'quasar_icon.png'}")

    icon_dark = render_icon_png(size=1024, on_dark=True)
    icon_dark.save(OUT / "quasar_icon_dark.png", "PNG")
    print(f"[png]   {OUT / 'quasar_icon_dark.png'}")

    # Full logos
    logo = render_full_logo(on_dark=False)
    logo.save(OUT / "quasar_logo.png", "PNG")
    print(f"[png]   {OUT / 'quasar_logo.png'}")

    logo_dark = render_full_logo(on_dark=True)
    logo_dark.save(OUT / "quasar_logo_dark.png", "PNG")
    print(f"[png]   {OUT / 'quasar_logo_dark.png'}")

    # Favicon: multi-resolution ICO
    sizes = [16, 32, 48, 64]
    favicon_imgs = [render_icon_png(size=s, on_dark=False) for s in sizes]
    favicon_imgs[0].save(OUT / "quasar_favicon.ico", format="ICO",
                          sizes=[(s, s) for s in sizes],
                          append_images=favicon_imgs[1:])
    print(f"[ico]   {OUT / 'quasar_favicon.ico'}")

    # Also save a wordmark-included SVG matching the PNG composition (for
    # consumers that prefer SVG over PNG).
    full_svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1600 480" width="1600" height="480">
  <defs>
    {ICON_SVG.split("<defs>")[1].split("</defs>")[0]}
  </defs>
  <rect width="1600" height="480" fill="#ffffff"/>
  <g transform="translate(86, 53) scale(0.37)">
    <use xlink:href="#__icon__"/>
  </g>
  <g transform="translate(86, 53)">
    <!-- inline icon (gradients reference defs above) -->
    <circle cx="190" cy="190" r="104" fill="url(#halo)"/>
    <polygon points="159,15 220,15 201,190 178,190" fill="url(#jet)"/>
    <polygon points="159,365 220,365 201,190 178,190" fill="url(#jet)"/>
    <ellipse cx="190" cy="190" rx="141" ry="24" fill="none" stroke="url(#disk)" stroke-width="5"/>
    <ellipse cx="190" cy="190" rx="96" ry="14" fill="none" stroke="url(#disk)" stroke-width="2" opacity="0.7"/>
    <circle cx="190" cy="190" r="52" fill="url(#core)"/>
    <circle cx="190" cy="190" r="8" fill="#ffffff"/>
  </g>
  <text x="540" y="245" font-family="Segoe UI, Inter, system-ui, sans-serif"
        font-size="180" font-weight="700" fill="#0b1e3a">Quasar</text>
  <text x="544" y="320" font-family="Segoe UI, Inter, system-ui, sans-serif"
        font-size="44" font-style="italic" fill="#475569">lymphoma structural-variant and fusion detection</text>
</svg>
'''
    write_svg(OUT / "quasar_logo.svg", full_svg)


if __name__ == "__main__":
    main()
