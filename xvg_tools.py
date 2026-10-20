"""GROMACS .xvg parsing and trajectory-metric summaries (stdlib only)."""
from __future__ import annotations

import statistics
from typing import Any


def read_xvg(text: str) -> dict[str, Any]:
    """Parse xvg text: skip '#' comments and '@' directives; extract title,
    axis labels, per-series legends, and numeric columns."""
    title = xlab = ylab = None
    legends: dict[int, str] = {}
    rows: list[list[float]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            body = line[1:].strip()
            if body.startswith("title") and len(body.split(None, 1)) > 1:
                title = body.split(None, 1)[1].strip().strip('"')
            elif body.startswith("xaxis") and len(body.split(None, 2)) > 2:
                xlab = body.split(None, 2)[2].strip('"')
            elif body.startswith("yaxis") and len(body.split(None, 2)) > 2:
                ylab = body.split(None, 2)[2].strip('"')
            elif body.startswith("s") and " legend " in body:
                try:
                    idx = int(body.split()[0][1:])
                    legends[idx] = body.split("legend", 1)[1].strip().strip('"')
                except (ValueError, IndexError):
                    pass
            continue
        parts = line.split()
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue  # skip malformed data lines
    if not rows:
        raise ValueError("no numeric data rows found in xvg content")
    width = min(len(r) for r in rows)
    return {
        "title": title,
        "x_label": xlab,
        "y_label": ylab,
        "series": [legends.get(i, f"col{i}") for i in range(width)],
        "n_points": len(rows),
        "columns": [[r[i] for r in rows] for i in range(width)],
    }


def _col_stats(v: list[float]) -> dict[str, Any]:
    return {
        "n": len(v),
        "first": v[0],
        "last": v[-1],
        "min": round(min(v), 5),
        "max": round(max(v), 5),
        "mean": round(statistics.fmean(v), 5),
        "stdev": round(statistics.stdev(v), 5) if len(v) > 1 else 0.0,
    }


def summary(text: str) -> dict[str, Any]:
    """Per-column descriptive stats for any .xvg file."""
    d = read_xvg(text)
    return {
        "title": d["title"],
        "x_label": d["x_label"],
        "y_label": d["y_label"],
        "n_points": d["n_points"],
        "columns": {d["series"][i]: _col_stats(col)
                    for i, col in enumerate(d["columns"])},
    }


def _mean_over(v: list[float], start: int, end: int) -> float:
    return statistics.fmean(v[start:end])


def rmsd_summary(text: str) -> dict[str, Any]:
    """Equilibration assessment from an RMSD .xvg (col 0 = time, col 1 = value)."""
    d = read_xvg(text)
    if len(d["columns"]) < 2:
        raise ValueError("expected >= 2 columns (time, RMSD) in the xvg file")
    t, y = d["columns"][0], d["columns"][1]
    n = len(y)
    k = max(1, n // 10)
    early = statistics.fmean(y[:k])
    late = statistics.fmean(y[-k:])
    drift_pct = (late - early) / max(abs(early), 1e-12) * 100.0
    tol = 0.1 * max(abs(late), 1e-12)
    plateau_idx = None
    for i in range(n - k):
        if (abs(y[i] - late) <= tol
                and all(abs(v - late) <= 1.5 * tol for v in y[i:i + k])):
            plateau_idx = i
            break
    return {
        "title": d["title"],
        "n_points": n,
        "time_range": [t[0], t[-1]],
        "value_stats": _col_stats(y),
        "early_mean_first10pct": round(early, 5),
        "late_mean_last10pct": round(late, 5),
        "drift_pct": round(drift_pct, 2),
        "plateau_reached": plateau_idx is not None,
        "plateau_time": round(t[plateau_idx], 5) if plateau_idx is not None else None,
        "verdict": ("appears equilibrated (drift < 10% between first and last 10%)"
                    if abs(drift_pct) < 10.0 else
                    "not clearly equilibrated (drift >= 10%) - inspect the trace"),
    }


def rmsf_peaks(text: str, top_n: int = 10) -> dict[str, Any]:
    """Most flexible residues/atoms from an RMSF .xvg (col 0 = index,
    col 1 = RMSF), plus contiguous flexible segments above mean + 1 stdev."""
    d = read_xvg(text)
    if len(d["columns"]) < 2:
        raise ValueError("expected >= 2 columns (index, RMSF) in the xvg file")
    x, y = d["columns"][0], d["columns"][1]
    n = len(y)
    mean = statistics.fmean(y)
    sd = statistics.stdev(y) if n > 1 else 0.0
    thresh = mean + sd
    order = sorted(range(n), key=lambda i: -y[i])[:max(1, min(top_n, n))]
    peaks = [{"index": x[i], "rmsf": round(y[i], 4)} for i in sorted(order)]
    segs, start, prev = [], None, None
    for i in range(n):
        if y[i] >= thresh:
            if start is None:
                start = i
            prev = i
        elif start is not None:
            segs.append((start, prev))
            start = None
    if start is not None:
        segs.append((start, prev))
    return {
        "title": d["title"],
        "n_points": n,
        "mean_rmsf": round(mean, 5),
        "threshold_used": round(thresh, 5),
        "most_flexible": peaks,
        "flexible_segments": [{"from_index": x[a], "to_index": x[b],
                               "n_points": b - a + 1} for a, b in segs[:20]],
    }
