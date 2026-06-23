# Quasar — Brand Assets

Logo & icon assets for **Quasar**, an NGS fusion-detection app for lymphoma.
The mark is a quasar: a luminous core, a doppler-bright accretion disk, and relativistic jets — a bright signal pulled from cosmic noise.

![Quasar](png/quasar-banner-dark.png)

## Files

### Vector (preferred — scales infinitely)
| File | Use |
|------|-----|
| `quasar-mark.svg` | Symbol only, tuned for **dark** backgrounds (transparent) |
| `quasar-mark-light.svg` | Symbol only, tuned for **light** backgrounds (transparent) |
| `quasar-icon.svg` | Full **app icon** — squircle deep-space tile, 512px |
| `quasar-mark-animated.svg` | **Animated** mark — pulsing core + orbiting accretion knots (great as a loading spinner) |
| `favicon.ico` | Multi-resolution favicon (16 / 32 / 48 px in one file) |

### Raster (`png/`)
| File | Use |
|------|-----|
| `quasar-icon-{16,32,48,64,128,180,192,256,512,1024}.png` | App icons / favicons / store listings |
| `quasar-mark-{256,512,1024}.png` | Symbol on dark (transparent) |
| `quasar-mark-light-{256,512}.png` | Symbol on light (transparent) |
| `quasar-banner-dark.png` / `quasar-banner-light.png` | README / social headers |
| `quasar-github-dark.png` / `quasar-github-light.png` | **GitHub headers** — backgrounds matched to GitHub's canvas (`#0d1117` / `#ffffff`) so they blend seamlessly |

## Common references

```html
<!-- Favicon -->
<link rel="icon" href="brand/favicon.ico" sizes="any">
<link rel="icon" type="image/png" sizes="32x32" href="brand/png/quasar-icon-32.png">
<link rel="apple-touch-icon" sizes="180x180" href="brand/png/quasar-icon-180.png">

<!-- Animated mark / loading spinner -->
<img src="brand/quasar-mark-animated.svg" width="120" alt="loading">
```

### Theme-aware GitHub header (auto dark/light)
Drop this at the top of your `README.md` — GitHub swaps the image to match the viewer's theme and it blends into the page background:

```html
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="brand/png/quasar-github-dark.png">
  <img alt="Quasar" src="brand/png/quasar-github-light.png" width="640">
</picture>

<!-- README header (GitHub) -->
<img src="brand/png/quasar-banner-dark.png" alt="Quasar" width="640">
```

## Brand tokens

| Token | Hex | Role |
|-------|-----|------|
| Cosmos | `#14162E` | Deep background / text |
| Disk Violet | `#6D54D1` | Disk outer / accent |
| Jet Cyan | `#1DA7C9` | Core / jet (light bg) |
| Bright Cyan | `#41C9EC` | Core / disk highlight (dark bg) |
| Slate | `#6B6E88` | Muted text / labels |

**Type** — Wordmark & headings: **Sora** (600). Labels & data: **IBM Plex Mono**.

> **Note:** `quasar-mark-animated.svg` animates when opened directly, embedded in an app, or used via `<img>`/`<object>`. GitHub strips SVG animation in rendered READMEs — use a PNG/GIF there if you need motion.

## Clearspace & minimums
- Keep clearspace ≥ the core-glow height on all sides.
- Minimum icon size: 24px. Minimum wordmark cap-height: 16px.
- Don't recolor the core, rotate the mark, or stretch the disk.
