"""GROMACS run-directory reports and cross-run comparison (stdlib only).

md_run_report replaces the one-off summarize/inspect scripts: one call over
a run directory gives per-xvg verdicts, crash detection, and a GO/NO-GO
style verdict, optionally recorded in the run registry.
"""
from __future__ import annotations

import os
from typing import Any

import run_registry
import xvg_tools

CRASH_MARKERS = ("fatal error", "segmentation fault", "sigsegv",
                 "water molecule starting", "not converging")


def _read(path: str) -> str:
    return open(path, errors="replace").read()


def _rmsd_in_dir(run_dir: str):
    """First parsable *rmsd*.xvg summary in a run directory, or None."""
    for name in sorted(os.listdir(run_dir)):
        if name.endswith(".xvg") and "rmsd" in name.lower():
            try:
                return xvg_tools.rmsd_summary(_read(os.path.join(run_dir, name)))
            except ValueError:
                continue
    return None


def md_run_report(run_dir: str, log: bool = True) -> dict[str, Any]:
    """Analyze a GROMACS run directory: xvg verdicts, log tail crash scan,
    artifact presence, and an overall verdict (GO / FINISHED / INCOMPLETE /
    CRASHED). Optionally records itself in the run registry."""
    if not os.path.isdir(run_dir):
        raise FileNotFoundError(f"run directory not found: {run_dir}")
    names = sorted(os.listdir(run_dir))
    xvgs = [n for n in names if n.endswith(".xvg")]

    rmsd_reports: list[dict] = []
    other_xvg: list[dict] = []
    for name in xvgs[:20]:
        path = os.path.join(run_dir, name)
        try:
            if "rmsd" in name.lower():
                rmsd_reports.append({"file": name,
                                     **xvg_tools.rmsd_summary(_read(path))})
            else:
                s = xvg_tools.summary(_read(path))
                other_xvg.append({"file": name, "n_points": s["n_points"],
                                  "columns": list(s["columns"])})
        except ValueError:
            continue

    log_file = next((n for n in names if n.endswith(".log")), None)
    crash = None
    finished = False
    if log_file:
        tail = _read(os.path.join(run_dir, log_file))[-20000:].lower()
        for marker in CRASH_MARKERS:
            if marker in tail:
                crash = marker
                break
        finished = "finished mdrun" in tail

    has_final = any(n in ("md.gro", "md_final.pdb") or n.endswith("_final.pdb")
                    for n in names)
    has_cpt = any(n.endswith(".cpt") for n in names)

    if crash:
        verdict = f"CRASHED (log marker: '{crash}')"
    elif rmsd_reports and any(r["plateau_reached"] for r in rmsd_reports) and finished:
        verdict = "GO - finished and RMSD plateau detected"
    elif finished:
        verdict = "FINISHED - but no clear RMSD plateau; inspect before using"
    else:
        verdict = "INCOMPLETE - no finished mdrun marker found"

    out = {"run_dir": run_dir, "verdict": verdict,
           "finished_mdrun": finished, "crash_marker": crash,
           "has_final_structure": has_final, "has_checkpoint": has_cpt,
           "n_xvg": len(xvgs), "rmsd_reports": rmsd_reports[:3],
           "other_xvg": other_xvg[:10]}
    if log:
        try:
            res = run_registry.log_run(
                kind="md_report",
                run_id=os.path.basename(os.path.normpath(run_dir)),
                status=verdict.split(" - ")[0],
                metrics={"verdict": verdict, "n_xvg": len(xvgs),
                         "crash": crash, "finished": finished},
                paths=[run_dir])
            out["registry_id"] = res["id"]
        except OSError:
            out["registry_id"] = None
    return out


def compare_runs(dir_a: str, dir_b: str) -> dict[str, Any]:
    """Side-by-side RMSD verdicts for two run directories (e.g. apo vs holo)."""
    ra, rb = _rmsd_in_dir(dir_a), _rmsd_in_dir(dir_b)
    if ra is None or rb is None:
        raise ValueError("could not find a parsable *rmsd*.xvg in one/both dirs")
    delta_late = ((rb["late_mean_last10pct"] - ra["late_mean_last10pct"])
                  / max(abs(ra["late_mean_last10pct"]), 1e-12) * 100.0)
    return {
        "dir_a": dir_a, "dir_b": dir_b,
        "a": {"plateau_reached": ra["plateau_reached"],
              "late_mean_last10pct": ra["late_mean_last10pct"],
              "drift_pct": ra["drift_pct"], "verdict": ra["verdict"]},
        "b": {"plateau_reached": rb["plateau_reached"],
              "late_mean_last10pct": rb["late_mean_last10pct"],
              "drift_pct": rb["drift_pct"], "verdict": rb["verdict"]},
        "late_mean_rel_difference_pct": round(delta_late, 2),
        "note": "late_mean_rel_difference_pct = (b - a)/a * 100 of the "
                "last-10% RMSD means",
    }
