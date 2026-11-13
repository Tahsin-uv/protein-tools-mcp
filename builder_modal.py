"""Modal-based detached, chunked bootstrap for the PDB metal-site index.

Implements BOOTSTRAP_SPEC.md: sync_batch / build_batch / merge_and_finalize
plus an orchestrator with pilot mode. Checkpoint = volume state (idempotent
by layout); failures are per-entry, never per-run.

Usage:
  modal run builder_modal.py --stage status
  modal run builder_modal.py --stage pilot            # 500 entries e2e
  modal run builder_modal.py --stage sync|build|merge|all
  modal run --detach builder_modal.py --stage all     # leave it running
"""
from __future__ import annotations

import json
import os

import modal

QC_VERSION = "v1"
SYNC_BATCH = 500
BUILD_BATCH = 250
MAX_MB = 50.0
LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))
LIST_JSON = os.path.join(LOCAL_DIR, "data", "pdb_index", "metal_entries.json")

img = modal.Image.debian_slim(python_version="3.11").add_local_python_source(
    "cif_tools", "pdb_utils", "builder_metal_index", "builder_modal")

app = modal.App("metal-index-bootstrap", image=img)
vol = modal.Volume.from_name("metal-index-vol", create_if_missing=True)

VOL = "/vol"
SUBSET = f"{VOL}/subset"
SHARDS = f"{VOL}/shards/build_{QC_VERSION}"
REPORTS = f"{VOL}/reports"
MERGED = f"{VOL}/merged"


def _commit():
    """Persist volume writes across containers (no-op on newer volumes)."""
    try:
        vol.commit()
    except Exception:
        pass


def _load_ids() -> list[str]:
    with open(LIST_JSON) as f:
        mapping = json.load(f)["mapping"]
    return sorted({i for v in mapping.values() for i in v})


def _batches(ids: list[str], size: int) -> list[tuple[str, list[str]]]:
    return [(f"{i:08d}", ids[i:i + size]) for i in range(0, len(ids), size)]


@app.function(image=img, volumes={VOL: vol}, timeout=7200, retries=3,
              max_containers=20)
def sync_batch(batch_id: str, ids: list[str]) -> dict:
    """Download one batch of entries into the volume (16 threads, resumable)."""
    from concurrent.futures import ThreadPoolExecutor
    import time as _t
    import urllib.error
    import urllib.request
    os.makedirs(SUBSET, exist_ok=True)

    def dl(eid: str) -> dict:
        dest = f"{SUBSET}/{eid}.cif.gz"
        if os.path.exists(dest):
            return {"id": eid, "result": "exists"}
        last: Exception | None = None
        for attempt in range(5):
            try:
                req = urllib.request.Request(
                    f"https://files.rcsb.org/download/{eid}.cif.gz",
                    headers={"User-Agent": "protein-tools-modal/0.1"})
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                if len(data) > MAX_MB * 1_000_000:
                    return {"id": eid, "result": "too_large",
                            "size_mb": round(len(data) / 1e6, 1)}
                with open(dest, "wb") as f:
                    f.write(data)
                return {"id": eid, "result": "ok",
                        "size_mb": round(len(data) / 1e6, 2)}
            except Exception as exc:
                last = exc
                _t.sleep(5 * (attempt + 1))
        return {"id": eid, "result": "error", "error": str(last)[:120]}

    with ThreadPoolExecutor(max_workers=16) as ex:
        results = list(ex.map(dl, ids))
    counts: dict[str, int] = {}
    for r in results:
        counts[r["result"]] = counts.get(r["result"], 0) + 1
    os.makedirs(REPORTS, exist_ok=True)
    with open(f"{REPORTS}/sync_{batch_id}.json", "w") as f:
        json.dump({"batch_id": batch_id, "n": len(ids), "counts": counts,
                   "results": results}, f)
    _commit()
    return {"stage": "sync", "batch_id": batch_id, "n": len(ids), **counts}
@app.function(image=img, volumes={VOL: vol}, timeout=7200, retries=3,
              memory=8192, max_containers=20)
def build_batch(batch_id: str, ids: list[str]) -> dict:
    """QC one batch of downloaded entries into a shard (sequential per
    container; scale comes from horizontal fan-out of many containers)."""
    import builder_metal_index as bmi
    os.makedirs(SHARDS, exist_ok=True)
    shard = f"{SHARDS}/batch_{batch_id}.jsonl"
    ok = shard + ".ok"
    if os.path.exists(ok):
        return {"stage": "build", "batch_id": batch_id, "result": "exists"}
    metals = sorted(set(bmi.DEFAULT_METALS))
    rows, errors, skipped = [], [], 0
    keys = [c.strip() for c in bmi.SCHEMA_COLS.split(",")]
    for eid in ids:
        path = f"{SUBSET}/{eid}.cif.gz"
        if not os.path.exists(path):
            errors.append({"entry_id": eid,
                           "error": "file missing (sync failed/skipped)"})
            continue
        try:
            _eid, site_rows, err = bmi._process_file(path, metals)
            if err:
                errors.append({"entry_id": eid, "error": err})
                continue
            for r in site_rows:
                d = dict(zip(keys, r))
                d["qc_version"] = QC_VERSION
                d["batch_id"] = batch_id
                rows.append(d)
        except Exception as exc:
            errors.append({"entry_id": eid, "error": f"{type(exc).__name__}: {exc}"[:150]})
    with open(shard, "w") as f:
        for d in rows:
            f.write(json.dumps({"type": "site", **d}, default=str) + "\n")
        for e in errors:
            f.write(json.dumps({"type": "error", **e}) + "\n")
    with open(ok, "w") as f:
        f.write(f"{len(rows)} sites, {len(errors)} errors, "
                f"{len(ids) - len(rows) - len(errors)} missing")
    _commit()
    return {"stage": "build", "batch_id": batch_id, "sites": len(rows),
            "errors": len(errors), "missing_files": skipped}


@app.function(image=img, volumes={VOL: vol}, timeout=7200, memory=8192)
def merge_and_finalize() -> dict:
    """Single-writer merge of all committed shards into metal_sites.db,
    followed by the validation queries from the spec."""
    import sqlite3
    os.makedirs(MERGED, exist_ok=True)
    db_path = f"{MERGED}/metal_sites.db"
    if os.path.exists(db_path):
        os.remove(db_path)  # merge is idempotent: always rebuild from shards
    shard_files = sorted(
        f for f in os.listdir(SHARDS) if f.endswith(".jsonl")
        and os.path.exists(os.path.join(SHARDS, f + ".ok")))
    con = sqlite3.connect(db_path)
    con.execute(
        "CREATE TABLE IF NOT EXISTS sites ("
        "entry_id TEXT, metal TEXT, site TEXT, coordination_number INTEGER, "
        "geometry TEXT, rmsd_deg REAL, n_imidazole INTEGER, n_carboxylate "
        "INTEGER, n_water INTEGER, n_backbone INTEGER, n_thiolate INTEGER, "
        "n_donors INTEGER, n_adjacent_metal INTEGER, bidentate INTEGER, "
        "artifact INTEGER, ligands TEXT, qc_version TEXT, batch_id TEXT)")
    n_rows = n_errors = n_shards = 0
    for sf in shard_files:
        with open(os.path.join(SHARDS, sf), errors="replace") as f:
            for line in f:
                rec = json.loads(line)
                if rec.get("type") == "site":
                    artifact = 1 if (rec.get("coordination_number") or 0) == 0 else 0
                    con.execute(
                        "INSERT INTO sites VALUES (?,?,?,?,?,?,?,?,?,?,?,"
                        "?,?,?,?,?,?,?)",
                        (rec["entry_id"], rec["metal"], rec["site"],
                         rec["coordination_number"], rec["geometry"],
                         rec["rmsd_deg"], rec["n_imidazole"],
                         rec["n_carboxylate"], rec["n_water"],
                         rec["n_backbone"], rec["n_thiolate"],
                         rec["n_donors"], rec["n_adjacent_metal"],
                         rec["bidentate"], artifact, rec["ligands"],
                         rec.get("qc_version", QC_VERSION),
                         rec.get("batch_id", "")))
                    n_rows += 1
                else:
                    n_errors += 1
        n_shards += 1
    con.execute("CREATE INDEX IF NOT EXISTS idx_metal ON sites(metal)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_geo ON sites(metal, geometry)")
    con.commit()
    stats = {
        "n_shards": n_shards, "n_site_rows": n_rows,
        "n_error_rows": n_errors,
        "n_entries": con.execute(
            "SELECT COUNT(DISTINCT entry_id) FROM sites").fetchone()[0],
        "per_metal": dict(con.execute(
            "SELECT metal, COUNT(*) FROM sites GROUP BY metal").fetchall()),
        "geometry_distribution": dict(con.execute(
            "SELECT geometry, COUNT(*) FROM sites GROUP BY geometry "
            "ORDER BY COUNT(*) DESC").fetchall()),
        "n_artifact_cn0": con.execute(
            "SELECT COUNT(*) FROM sites WHERE artifact = 1").fetchone()[0],
        "n_cn1_unusual": con.execute(
            "SELECT COUNT(*) FROM sites WHERE coordination_number = 1"
        ).fetchone()[0],
    }
    con.close()
    report = {"validated": True,  # flag-don't-fail: artifacts are classified,
              # not fatal - query tools default to artifact=0 AND cn >= 2
              "classification_note": ("CN-0 rows flagged artifact=1 "
                                      "(isolated-ion artifact); CN-1 rows "
                                      "kept as unusual chemistry"),
              "db_path": db_path, **stats}
    os.makedirs(REPORTS, exist_ok=True)
    with open(f"{REPORTS}/merge_final.json", "w") as f:
        json.dump(report, f, indent=1)
    _commit()
    return report


@app.function(image=img, volumes={VOL: vol}, timeout=600)
def volume_inventory() -> dict:
    """Walk the volume and report per-directory file counts and byte sizes."""
    import collections
    inv: dict[str, dict] = {}
    total_files = 0
    total_bytes = 0
    for root, _dirs, files in os.walk(VOL):
        rel = os.path.relpath(root, VOL)
        key = rel.split(os.sep)[0]
        n = len(files)
        b = 0
        exts: dict[str, int] = collections.Counter()
        for f in files:
            fp = os.path.join(root, f)
            try:
                b += os.path.getsize(fp)
            except OSError:
                pass
            ext = os.path.splitext(f)[1] or "(noext)"
            exts[ext] += 1
        if key not in inv:
            inv[key] = {"files": 0, "bytes": 0, "exts": {}}
        inv[key]["files"] += n
        inv[key]["bytes"] += b
        for k, v in exts.items():
            inv[key]["exts"][k] = inv[key]["exts"].get(k, 0) + v
        total_files += n
        total_bytes += b
    return {"total_files": total_files, "total_gib": round(total_bytes / (1024 ** 3), 2),
            "dirs": inv}


@app.function(image=img, volumes={VOL: vol}, timeout=600)
def status() -> dict:
    def count_dir(d: str, *exts: str) -> int:
        if not os.path.isdir(d):
            return 0
        return len([f for f in os.listdir(d)
                    if not exts or any(f.endswith(e) for e in exts)])
    return {"qc_version": QC_VERSION,
            "synced_files": count_dir(SUBSET, ".cif.gz"),
            "shards": count_dir(SHARDS, ".jsonl"),
            "shards_committed": count_dir(SHARDS, ".jsonl.ok"),
            "reports": count_dir(REPORTS),
            "db_exists": os.path.exists(f"{MERGED}/metal_sites.db")}
@app.local_entrypoint()
def main(stage: str = "pilot", limit: int = SYNC_BATCH):
    """Orchestrator.

    stage: pilot (first `limit` entries end-to-end) | sync | build | merge
           | all | status | inventory
    """
    ids = _load_ids()
    print(f"entry map: {len(ids)} unique entries")
    if stage == "status":
        print(json.dumps(status.remote(), indent=1))
        return
    if stage == "inventory":
        print(json.dumps(volume_inventory.remote(), indent=1))
        return
    if stage not in ("pilot", "sync", "build", "merge", "all"):
        raise SystemExit(f"unknown stage: {stage}")
    if stage == "pilot":
        ids = ids[:limit]
        print(f"PILOT: running end-to-end on {len(ids)} entries")
    if stage in ("pilot", "sync", "all"):
        batches = _batches(ids, SYNC_BATCH)
        print(f"[sync] {len(batches)} batches x {SYNC_BATCH} ids")
        for r in sync_batch.starmap(batches):
            print(r)
    if stage in ("pilot", "build", "all"):
        batches = _batches(ids, BUILD_BATCH)
        print(f"[build] {len(batches)} batches x {BUILD_BATCH} ids")
        for r in build_batch.starmap(batches):
            print(r)
    if stage in ("pilot", "merge", "all"):
        print("[merge]")
        print(json.dumps(merge_and_finalize.remote(), indent=1))
    print("done")


