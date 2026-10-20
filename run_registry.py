"""SQLite-backed run registry: cross-session memory for MD/project runs.

This is the piece ad-hoc scripts cannot provide: a persistent, queryable
record of what was run, when, with what status and metrics.
"""
from __future__ import annotations

import datetime
import json
import os
import sqlite3
from typing import Any

import app_paths

_DB = str(app_paths.state_file("runs", "runs.db")[0])


def db_path() -> str:
    """Filesystem path of the registry DB (per-user app-data dir)."""
    return _DB


def _conn() -> sqlite3.Connection:
    os.makedirs(os.path.dirname(_DB), exist_ok=True)
    con = sqlite3.connect(_DB)
    con.execute(
        "CREATE TABLE IF NOT EXISTS runs ("
        "id INTEGER PRIMARY KEY AUTOINCREMENT, "
        "ts TEXT NOT NULL, "
        "run_id TEXT, "
        "kind TEXT NOT NULL, "
        "status TEXT, "
        "notes TEXT, "
        "metrics TEXT, "
        "paths TEXT)")
    return con


def _dump(x) -> str:
    return json.dumps(x if x is not None else {}, default=str)


def log_run(kind: str, run_id: str | None = None, status: str = "logged",
            notes: str = "", metrics: Any = None,
            paths: Any = None) -> dict[str, Any]:
    """Insert a run/event record; returns its registry id and timestamp."""
    con = _conn()
    ts = datetime.datetime.now().isoformat(timespec="seconds")
    cur = con.execute(
        "INSERT INTO runs (ts, run_id, kind, status, notes, metrics, paths) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (ts, run_id, kind, status, notes, _dump(metrics), _dump(paths)))
    con.commit()
    rid = cur.lastrowid
    con.close()
    return {"id": rid, "ts": ts, "kind": kind, "run_id": run_id,
            "status": status}


def list_runs(kind: str | None = None, limit: int = 20) -> dict[str, Any]:
    """Most recent registry entries, newest first, optionally by kind."""
    con = _conn()
    q = ("SELECT id, ts, run_id, kind, status FROM runs "
         + ("WHERE kind = ? " if kind else "")
         + "ORDER BY id DESC LIMIT ?")
    args = (kind, limit) if kind else (limit,)
    rows = con.execute(q, args).fetchall()
    con.close()
    return {"db": _DB,
            "n": len(rows),
            "runs": [{"id": r[0], "ts": r[1], "run_id": r[2], "kind": r[3],
                      "status": r[4]} for r in rows]}


def run_details(rid: str) -> dict[str, Any]:
    """Full record for a registry entry, by numeric id or run_id (latest)."""
    con = _conn()
    if str(rid).isdigit():
        row = con.execute("SELECT * FROM runs WHERE id = ?", (int(rid),)
                          ).fetchone()
    else:
        row = con.execute(
            "SELECT * FROM runs WHERE run_id = ? ORDER BY id DESC LIMIT 1",
            (str(rid),)).fetchone()
    con.close()
    if row is None:
        raise ValueError(f"no registry entry for '{rid}'")
    cols = ["id", "ts", "run_id", "kind", "status", "notes", "metrics",
            "paths"]
    d = dict(zip(cols, row))
    for k in ("metrics", "paths"):
        try:
            d[k] = json.loads(d[k]) if d[k] else None
        except (ValueError, TypeError):
            pass
    return d


def update_run(rid: int, status: str | None = None, notes: str | None = None,
               metrics: Any = None) -> dict[str, Any]:
    """Update status/notes/metrics of an existing entry."""
    con = _conn()
    sets, args = [], []
    if status is not None:
        sets.append("status = ?")
        args.append(status)
    if notes is not None:
        sets.append("notes = ?")
        args.append(notes)
    if metrics is not None:
        sets.append("metrics = ?")
        args.append(_dump(metrics))
    if not sets:
        con.close()
        return {"updated": 0}
    args.append(int(rid))
    con.execute(f"UPDATE runs SET {', '.join(sets)} WHERE id = ?", args)
    con.commit()
    n = con.total_changes
    con.close()
    return {"updated": n, "id": int(rid)}
