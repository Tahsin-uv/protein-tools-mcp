"""Plotting tools for GROMACS .xvg data (matplotlib, Agg backend).

matplotlib is imported LAZILY inside the plot functions so the MCP server
starts (and all non-plotting tools work) even on machines where matplotlib
is missing or its native DLL is blocked (e.g. Windows Application Control).
"""
from __future__ import annotations

import os
from typing import Any

import xvg_tools

MAX_POINTS = 5000  # decimate larger series before plotting

_PLT_OK: bool | None = None


def _matplotlib_usable() -> bool:
    """Probe matplotlib ONCE in a subprocess with a hard timeout - on some
    locked-down machines the import can hang (e.g. DLL policy), which must
    never hang the whole MCP server."""
    global _PLT_OK
    if _PLT_OK is None:
        import subprocess
        import sys
        try:
            r = subprocess.run(
                [sys.executable, "-c",
                 "import matplotlib; matplotlib.use('Agg'); "
                 "import matplotlib.pyplot"],
                capture_output=True, timeout=15)
            _PLT_OK = r.returncode == 0
        except Exception:
            _PLT_OK = False
    return _PLT_OK


def _plt():
    """Lazy matplotlib import; raises a clear error when unavailable."""
    if not _matplotlib_usable():
        raise RuntimeError(
            "matplotlib is unavailable or blocked on this machine - install "
            "it with 'pip install matplotlib' to use the plotting tools")
    import matplotlib
    matplotlib.use("Agg")  # headless: no display server needed
    import matplotlib.pyplot as plt
    return plt


def _decimate(values: list, max_points: int = MAX_POINTS) -> list[int]:
    if len(values) <= max_points:
        return list(range(len(values)))
    step = len(values) // max_points + 1
    return list(range(0, len(values), step))


def _file_size(path: str) -> int:
    return os.path.getsize(path) if os.path.exists(path) else 0


def plot_xvg(xvg_path: str, out_png: str = "", title: str = "",
             y_from_zero: bool = False) -> dict[str, Any]:
    """Plot every series of an .xvg file against column 0 into a PNG."""
    d = xvg_tools.read_xvg(open(xvg_path, errors="replace").read())
    if not out_png:
        out_png = os.path.splitext(xvg_path)[0] + ".png"
    x = d["columns"][0]
    idx = _decimate(x)
    plt = _plt()
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=130)
    for i in range(1, len(d["columns"])):
        y = d["columns"][i]
        label = d["series"][i] if i < len(d["series"]) else f"col{i}"
        ax.plot([x[j] for j in idx], [y[j] for j in idx], lw=1.0, label=label)
    ax.set_xlabel(d["x_label"] or "x")
    ax.set_ylabel(d["y_label"] or "")
    ax.set_title(title or d["title"] or os.path.basename(xvg_path))
    if len(d["columns"]) > 2:
        ax.legend(fontsize=8)
    if y_from_zero:
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return {"out_png": out_png, "size_bytes": _file_size(out_png),
            "series": d["series"][1:], "points_plotted": len(idx),
            "n_points_total": d["n_points"]}


def compare_plot(xvg_paths: str, out_png: str = "", labels: str = "",
                 value_column: int = 1, title: str = "",
                 y_from_zero: bool = False) -> dict[str, Any]:
    """Overlay the value column of several .xvg files (comma-separated
    paths, e.g. apo vs holo RMSD) into one PNG for direct comparison."""
    paths = [p.strip() for p in xvg_paths.split(",") if p.strip()]
    if len(paths) < 2:
        raise ValueError("compare_plot needs at least two .xvg paths "
                         "(comma-separated)")
    if not out_png:
        out_png = os.path.join(os.path.dirname(paths[0]), "compare_plot.png")
    if labels:
        labs = [l.strip() for l in labels.split(",")]
    else:
        labs = [os.path.basename(p) for p in paths]
    plt = _plt()
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=130)
    plotted = []
    for i, p in enumerate(paths):
        d = xvg_tools.read_xvg(open(p, errors="replace").read())
        if value_column >= len(d["columns"]):
            raise ValueError(f"{p}: column {value_column} does not exist "
                             f"(has {len(d['columns'])} columns)")
        x, y = d["columns"][0], d["columns"][value_column]
        idx = _decimate(x)
        ax.plot([x[j] for j in idx], [y[j] for j in idx], lw=1.0,
                label=labs[i] if i < len(labs) else os.path.basename(p))
        plotted.append({"file": os.path.basename(p), "series": d["series"],
                        "n_points": d["n_points"]})
    ax.set_xlabel(d["x_label"] or "x")
    ax.set_ylabel(d["y_label"] or "")
    ax.set_title(title or "comparison")
    ax.legend(fontsize=8)
    if y_from_zero:
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return {"out_png": out_png, "size_bytes": _file_size(out_png),
            "plotted": plotted, "value_column": value_column}
