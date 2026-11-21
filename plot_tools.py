"""Plotting tools for GROMACS .xvg data (matplotlib, Agg backend).

matplotlib is imported LAZILY inside the plot functions so the MCP server
starts (and all non-plotting tools work) even on machines where matplotlib
is missing or its native DLL is blocked (e.g. Windows Application Control).
"""
from __future__ import annotations

import os
import re
import stat
from pathlib import Path
from typing import Any

import xvg_tools

MAX_POINTS = 5000  # decimate larger series before plotting

# Only PNG output is permitted — prevents writing unexpected file types.
PLOT_ALLOWED_EXTENSIONS = {".png"}

# Windows reserved device names (case-insensitive, with or without extension).
_WINDOWS_RESERVED = {
    "CON", "PRN", "AUX", "NUL",
    "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9",
    "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9",
}

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


def _plot_output_dir() -> Path:
    """Directory that plot PNGs may be written to.

    Configurable via ``PROTEIN_TOOLS_PLOT_DIR``; defaults to the current
    working directory when unset. Created on demand.
    """
    override = os.environ.get("PROTEIN_TOOLS_PLOT_DIR", "").strip()
    if override:
        d = Path(override).expanduser().resolve()
    else:
        d = Path.cwd().resolve()
    d.mkdir(parents=True, exist_ok=True)
    return d


def _has_symlink_component(path: Path) -> bool:
    """Return True if any component of *path* (other than the final
    component) is a symlink.  This prevents symlink attacks where an
    attacker swaps a directory for a symlink after resolution but before
    the file is written."""
    parts = path.parts
    # Check every ancestor (not the final component, which will be created).
    for i in range(1, len(parts)):
        ancestor = Path(*parts[:i])
        try:
            if ancestor.is_symlink():
                return True
        except OSError:
            return True
    return False


def _is_windows_reserved_name(path: Path) -> bool:
    """Return True if the final path component is a Windows reserved
    device name (e.g. CON, COM1, LPT1), which would silently create
    or open a device instead of a file."""
    stem = path.stem.upper()
    return stem in _WINDOWS_RESERVED


def _safe_plot_output_path(user_path: str, default_name: str) -> Path:
    """Resolve and validate a plot output path.

    - Blocks ``..`` path traversal.
    - Restricts output to ``.png`` only.
    - Confines the result to the configurable output directory
      (``PROTEIN_TOOLS_PLOT_DIR`` env var, default cwd).
    - Blocks symlink components in the path (prevents symlink swaps).
    - Blocks Windows reserved device names.
    - Creates parent directories within the allowed dir.
    """
    out_dir = _plot_output_dir()
    p = Path(user_path)
    if ".." in p.parts:
        raise ValueError(
            f"plot output path may not contain '..': {user_path!r}")
    # Extension must be .png (case-insensitive).
    if p.suffix.lower() not in PLOT_ALLOWED_EXTENSIONS:
        raise ValueError(
            f"plot output must be a .png file, got {p.suffix!r}: "
            f"{user_path!r}")
    if _is_windows_reserved_name(p):
        raise ValueError(
            f"plot output name is a Windows reserved device: "
            f"{user_path!r}")
    if not p.is_absolute():
        p = out_dir / p
        # Re-check suffix after joining with out_dir.
        if p.suffix.lower() not in PLOT_ALLOWED_EXTENSIONS:
            raise ValueError(
                f"plot output must be a .png file, got {p.suffix!r}: "
                f"{user_path!r}")
        if _is_windows_reserved_name(p):
            raise ValueError(
                f"plot output name is a Windows reserved device: "
                f"{user_path!r}")
        if _has_symlink_component(p):
            raise ValueError(
                f"plot output path contains a symlink component: "
                f"{user_path!r}")
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    resolved = p.resolve()
    try:
        resolved.relative_to(out_dir)
    except ValueError:
        raise ValueError(
            f"plot output path is outside the allowed directory "
            f"{out_dir}: {user_path!r}")
    if _has_symlink_component(resolved):
        raise ValueError(
            f"plot output path contains a symlink component: "
            f"{user_path!r}")
    if _is_windows_reserved_name(resolved):
        raise ValueError(
            f"plot output name is a Windows reserved device: "
            f"{user_path!r}")
    p = resolved
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _safe_plot_input_path(user_path: str) -> Path:
    """Validate an input .xvg path: block ``..`` traversal and resolve
    symlinks so we know exactly which file is being read."""
    p = Path(user_path)
    if ".." in p.parts:
        raise ValueError(
            f"plot input path may not contain '..': {user_path!r}")
    resolved = p.resolve()
    # After resolution, re-check that no ancestor of the resolved path is
    # a symlink (defense against symlink-swap races).
    if _has_symlink_component(resolved):
        raise ValueError(
            f"plot input path contains a symlink component: {user_path!r}")
    return resolved


def plot_xvg(xvg_path: str, out_png: str = "", title: str = "",
             y_from_zero: bool = False) -> dict[str, Any]:
    """Plot every series of an .xvg file against column 0 into a PNG."""
    xvg_path = str(_safe_plot_input_path(xvg_path))
    d = xvg_tools.read_xvg(open(xvg_path, errors="replace").read())
    if not out_png:
        out_png_path = _plot_output_dir() / (
            os.path.splitext(os.path.basename(xvg_path))[0] + ".png")
    else:
        out_png_path = _safe_plot_output_path(
            out_png, os.path.basename(xvg_path) + ".png")
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
    fig.savefig(out_png_path)
    plt.close(fig)
    return {"out_png": str(out_png_path), "size_bytes": _file_size(str(out_png_path)),
            "series": d["series"][1:], "points_plotted": len(idx),
            "n_points_total": d["n_points"]}


def compare_plot(xvg_paths: str, out_png: str = "", labels: str = "",
                 value_column: int = 1, title: str = "",
                 y_from_zero: bool = False) -> dict[str, Any]:
    """Overlay the value column of several .xvg files (comma-separated
    paths, e.g. apo vs holo RMSD) into one PNG for direct comparison."""
    paths = [str(_safe_plot_input_path(p.strip())) for p in xvg_paths.split(",") if p.strip()]
    if len(paths) < 2:
        raise ValueError("compare_plot needs at least two .xvg paths "
                         "(comma-separated)")
    if not out_png:
        out_png_path = _plot_output_dir() / "compare_plot.png"
    else:
        out_png_path = _safe_plot_output_path(out_png, "compare_plot.png")
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
    fig.savefig(out_png_path)
    plt.close(fig)
    return {"out_png": str(out_png_path), "size_bytes": _file_size(str(out_png_path)),
            "plotted": plotted, "value_column": value_column}
