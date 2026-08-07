#!/usr/bin/env python3
"""Shared figure-saving helper: write every graph as both a raster PNG and a vector SVG.

Central so every plotting script in the pipeline emits publication-ready vector figures alongside
the raster ones with a single call. Use in place of ``fig.savefig(path, ...)``:

    from plot_utils import save_figure
    save_figure(fig, "foo.png", dpi=300, bbox_inches="tight")   # writes foo.png AND foo.svg
"""
import os

# Formats every graph is written in. PNG first (also the one embedded/referenced by the report).
FORMATS = ("png", "svg")


def save_figure(fig, path, formats=FORMATS, **savefig_kwargs):
    """Save ``fig`` to ``path`` in each of ``formats`` (default PNG + SVG).

    ``path`` may carry any image extension (or none); the stem is reused for every format, so
    ``save_figure(fig, "x/plot.png")`` writes ``x/plot.png`` and ``x/plot.svg``. ``dpi`` is dropped
    for vector formats (it is meaningless there). Returns the list of written paths.
    """
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
