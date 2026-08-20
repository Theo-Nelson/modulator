#!/usr/bin/env python3
"""Shared figure styling + saving helpers.

Central so every plotting script in the pipeline renders with one house style (Arial / Arial-like
sans-serif, enlarged fonts) and emits publication-ready vector figures (PDF + SVG) alongside the
raster PNG with a single call. Use in place of ``fig.savefig(path, ...)``:

    from plot_utils import save_figure
    save_figure(fig, "foo.png", dpi=300, bbox_inches="tight")   # writes foo.png, foo.pdf, foo.svg

``save_figure`` applies the house style (font family + FONT_BUMP) to the figure before writing, so
both the persisted files and any subsequent inline embed of the same figure share the look.
"""
import os

# Formats every graph is written in. PNG first (the one embedded/referenced by the report); PDF for
# dropping the graph straight into other figures/manuscripts; SVG kept for fully-editable vectors.
FORMATS = ("png", "pdf", "svg")

# Every font size in a figure is raised by this many points (publication-scale, per request). The
# title may then need two lines -- bbox_inches="tight" grows the canvas to fit.
FONT_BUMP = 12

# Arial first, then metric-compatible / Arial-like fallbacks (Nimbus Sans is URW's Helvetica/Arial
# clone and ships on most Linux boxes); DejaVu Sans is matplotlib's always-present last resort.
_SANS = ["Arial", "Helvetica", "Nimbus Sans", "Liberation Sans", "Arimo", "DejaVu Sans"]

_STYLE_APPLIED = False


def setup_matplotlib_style():
    """Set the process-wide house style (Arial/Arial-like sans-serif). Idempotent + cheap."""
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return
    try:
        import matplotlib
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = _SANS + list(
            matplotlib.rcParams.get("font.sans-serif", [])
        )
        # keep mathtext/minus consistent with the sans face
        matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
        matplotlib.rcParams["axes.unicode_minus"] = False
        matplotlib.rcParams["svg.fonttype"] = "none"   # keep text as text in the SVG/PDF, not paths
        matplotlib.rcParams["pdf.fonttype"] = 42        # embed TrueType so text stays editable
        _STYLE_APPLIED = True
    except Exception:
        pass


def bump_fonts(fig, delta=FONT_BUMP):
    """Raise every text element's font size in ``fig`` by ``delta`` points (idempotent per figure).

    Applied centrally so scripts keep their relative sizing (title > label > annotation) and every
    size shifts up together, rather than editing dozens of hard-coded ``fontsize=`` call sites.
    """
    if delta == 0 or getattr(fig, "_fonts_bumped", False):
        return
    try:
        from matplotlib.text import Text
        for t in fig.findobj(Text):
            try:
                sz = t.get_fontsize()
                if sz:
                    t.set_fontsize(sz + delta)
            except Exception:
                continue
        fig._fonts_bumped = True
    except Exception:
        pass
    # Reflow with the enlarged fonts. constrained_layout runs at DRAW time (i.e. after this bump),
    # so it re-allocates space for the now-larger titles/labels/ticks -- unlike a build-time
    # tight_layout(), which was computed against the old sizes and would otherwise let elements
    # collide. Skip figures that already declared an engine (e.g. layout="constrained" at build).
    try:
        if fig.get_layout_engine() is None:
            fig.set_layout_engine("constrained")
    except Exception:
        pass


def save_figure(fig, path, formats=FORMATS, **savefig_kwargs):
    """Save ``fig`` to ``path`` in each of ``formats`` (default PNG + PDF + SVG).

    Applies the house style + FONT_BUMP first. ``path`` may carry any image extension (or none);
    the stem is reused for every format, so ``save_figure(fig, "x/plot.png")`` writes
    ``x/plot.png``, ``x/plot.pdf`` and ``x/plot.svg``. ``dpi`` is dropped for vector formats.
    Returns the list of written paths.
    """
    setup_matplotlib_style()
    bump_fonts(fig)
    stem, ext = os.path.splitext(path)
    if ext.lower() not in (".png", ".svg", ".pdf", ".jpg", ".jpeg", ".eps", ".tif", ".tiff", ""):
        stem = path  # unrecognized extension -> treat the whole thing as the stem
    out_dir = os.path.dirname(stem)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    written = []
    for fmt in formats:
        kw = dict(savefig_kwargs)
        if fmt in ("svg", "pdf", "eps"):
            kw.pop("dpi", None)
        p = f"{stem}.{fmt}"
        fig.savefig(p, format=fmt, **kw)
        written.append(p)
    return written
