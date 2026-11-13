"""PDB-wide metal-site index builder (three resumable stages).

Pipeline:
  list   - RCSB Search API -> data/pdb_index/metal_entries.json
           (metal -> entry-id map, paginated, experimental structures only)
  sync   - files.rcsb.org  -> data/pdb_subset/*.cif.gz
           (parallel, resumable, size-capped)
  build  - multiprocess QC -> data/index/metal_sites.db
           (runs the tested metal_geometry verdicts over every entry)
  status - inventory of what exists where
  query  - geometry queries over the index (e.g. octahedral ZN sites)

All stages are idempotent: re-running skips completed work, so the first
bootstrap is a few hours of mostly unattended transfer and every later
update is minutes. Designed for Windows (spawn-safe multiprocessing).
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sqlite3
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import cif_tools
import pdb_utils

BASE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE, "data")
LIST_JSON = os.path.join(DATA_DIR, "pdb_index", "metal_entries.json")
DOWNLOAD_DIR = os.path.join(DATA_DIR, "pdb_subset")
DB_PATH = os.path.join(DATA_DIR, "index", "metal_sites.db")

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{eid}.cif.gz"
DEFAULT_METALS = ["ZN", "FE", "MG", "MN", "CU", "NI", "CO", "CA", "CD", "HG",
                  "HEM", "SF4", "FES", "F3S", "CLA", "CLB", "B12", "CNC"]


def _post_search(metal: str, start: int, rows: int, retries: int = 3) -> dict:
    import time
    import urllib.error
    import urllib.request
    body = {
        "query": {"type": "terminal", "service": "text",
                  "parameters": {
                      "attribute": ("rcsb_nonpolymer_entity_container_"
                                    "identifiers.nonpolymer_comp_id"),
                      "operator": "exact_match",
                      "value": metal.upper()}},
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "results_content_type": ["experimental"]},
        "return_type": "entry",
    }
    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(
                SEARCH_URL, data=json.dumps(body).encode("utf-8"),
                headers={"Content-Type": "application/json",
                         "User-Agent": "protein-tools-index/0.1"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            if attempt < retries:
                time.sleep(5 * (attempt + 1))
    raise last


def fetch_metal_entry_list(metals: list[str], out_json: str = LIST_JSON,
                           limit_per_metal: int | None = None) -> dict:
    """Stage 1: metal -> sorted entry-id map via the RCSB Search API.

    Saves progress after each metal, so interrupted runs keep completed
    metals and re-running skips them."""
    mapping: dict[str, list] = {}
    if os.path.exists(out_json):  # resume: keep already-fetched metals
        try:
            with open(out_json) as f:
                mapping = json.load(f).get("mapping", {})
        except (ValueError, OSError):
            mapping = {}
    for m in metals:
        if mapping.get(m.upper()):
            print(f"  {m.upper()}: {len(mapping[m.upper()])} entries "
                  f"(already fetched, skipping)")
            continue
        ids: list[str] = []
        start, rows = 0, 10000
        while True:
            resp = _post_search(m, start, rows)
            total = resp.get("total_count", 0)
            got = [x["identifier"] for x in resp.get("result_set", [])]
            ids.extend(got)
            start += len(got)
            if not got or start >= total:
                break
            if limit_per_metal and len(ids) >= limit_per_metal:
                break
        ids = sorted(set(ids))
        if limit_per_metal:
            ids = ids[:limit_per_metal]
        mapping[m.upper()] = ids
        print(f"  {m.upper()}: {len(ids)} entries"
              + (f" (of {total} total)" if limit_per_metal else ""))
        # save after every metal: a later failure keeps finished work
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as f:
            json.dump({"metals": {k: len(v) for k, v in mapping.items()},
                       "mapping": mapping}, f)
    n_all = len({i for v in mapping.values() for i in v})
    print(f"saved {out_json} ({n_all} unique entries)")
    return mapping
def _load_mapping(list_json: str) -> dict:
    with open(list_json) as f:
        return json.load(f)["mapping"]


def _download_one(eid: str, download_dir: str, max_mb: float) -> dict:
    import urllib.request
    dest = os.path.join(download_dir, f"{eid}.cif.gz")
    if os.path.exists(dest):
        return {"id": eid, "result": "exists"}
    try:
        req = urllib.request.Request(
            DOWNLOAD_URL.format(eid=eid),
            headers={"User-Agent": "protein-tools-index/0.1"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        if len(data) > max_mb * 1_000_000:
            return {"id": eid, "result": "too_large",
                    "size_mb": round(len(data) / 1e6, 1)}
        with open(dest, "wb") as f:
            f.write(data)
        return {"id": eid, "result": "ok", "size_mb": round(len(data) / 1e6, 2)}
    except Exception as exc:  # per-entry failure must not kill the batch
        return {"id": eid, "result": "error", "error": str(exc)[:120]}


def sync_entries(list_json: str = LIST_JSON, download_dir: str = DOWNLOAD_DIR,
                 workers: int = 8, max_mb: float = 50.0,
                 limit: int | None = None) -> dict:
    """Stage 2: resumable parallel download of the mapped entries."""
    mapping = _load_mapping(list_json)
    ids = sorted({i for v in mapping.values() for i in v})
    if limit:
        ids = ids[:limit]
    os.makedirs(download_dir, exist_ok=True)
    print(f"sync: {len(ids)} entries, {workers} workers, "
          f"max {max_mb} MB per file")
    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for res in pool.map(lambda e: _download_one(e, download_dir, max_mb),
                            ids):
            results.append(res)
    counts = Counter(r["result"] for r in results)
    failed = [r for r in results if r["result"] == "error"]
    skipped = [r for r in results if r["result"] == "too_large"]
    out = {"total": len(ids), **{f"n_{k}": v for k, v in counts.items()},
           "failed": failed[:10], "skipped_too_large": skipped[:10]}
    with open(os.path.join(download_dir, "_sync_report.json"), "w") as f:
        json.dump({"results": results}, f)
    print(f"sync done: {dict(counts)}")
    return out
SCHEMA_COLS = ("entry_id, metal, site, coordination_number, geometry, "
               "rmsd_deg, n_imidazole, n_carboxylate, n_water, n_backbone, "
               "n_thiolate, n_donors, n_adjacent_metal, bidentate, artifact, "
               "ligands")
INSERT_SQL = (f"INSERT INTO sites ({SCHEMA_COLS}) VALUES "
              "(" + ",".join("?" * len(SCHEMA_COLS.split(","))) + ")")


def _init_db(db_path: str) -> sqlite3.Connection:
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    con = sqlite3.connect(db_path)
    con.execute(
        "CREATE TABLE IF NOT EXISTS sites ("
        "entry_id TEXT, metal TEXT, site TEXT, "
        "coordination_number INTEGER, geometry TEXT, rmsd_deg REAL, "
        "n_imidazole INTEGER, n_carboxylate INTEGER, n_water INTEGER, "
        "n_backbone INTEGER, n_thiolate INTEGER, n_donors INTEGER, "
        "n_adjacent_metal INTEGER, bidentate INTEGER, "
        "artifact INTEGER, ligands TEXT)")
    # migration for indexes built with the older (13-column) schema
    cols = {r[1] for r in con.execute("PRAGMA table_info(sites)")}
    if "n_adjacent_metal" not in cols:
        con.execute("ALTER TABLE sites ADD COLUMN n_adjacent_metal INTEGER")
    if "bidentate" not in cols:
        con.execute("ALTER TABLE sites ADD COLUMN bidentate INTEGER")
    if "artifact" not in cols:
        con.execute("ALTER TABLE sites ADD COLUMN artifact INTEGER")
    con.execute("CREATE INDEX IF NOT EXISTS idx_sites_metal "
                "ON sites(metal)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_sites_geo "
                "ON sites(metal, geometry)")
    con.commit()
    return con


def _process_file(path: str, metals: list[str]):
    """Worker: parse one downloaded CIF and return (entry_id, rows, error).

    Top-level function: required for Windows spawn multiprocessing."""
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", errors="replace") as f:
            text = f.read()
        atoms = cif_tools.to_pdb_atoms(cif_tools.parse_cif(text))
        entry_id = os.path.basename(path).replace(".cif.gz", "").replace(".cif", "")
        present = sorted({a["element"].strip().upper() for a in atoms["atoms"]}
                         & set(metals))
        rows = []
        for m in present:
            g = pdb_utils.metal_geometry(atoms, m)
            for site in g["metal_sites"]:
                fit = site.get("geometry_fit") or {}
                types = Counter(c.get("type", "other")
                                for c in site["ligand_contacts"])
                n_adj = sum(1 for c in site["ligand_contacts"]
                            if str(c.get("type", "")).startswith("adjacent metal"))
                ligands = ";".join(f"{c['residue']}:{c['atom']}"
                                   for c in site["ligand_contacts"])
                rows.append((
                    entry_id, m, site["site"], site["coordination_number"],
                    fit.get("best"), fit.get("rmsd_deg"),
                    types.get("imidazole N", 0), types.get("carboxylate O", 0),
                    types.get("water", 0), types.get("backbone carbonyl O", 0),
                    types.get("thiol/thiolate S", 0),
                    len(site["ligand_contacts"]), n_adj,
                    1 if site.get("bidentate_chelates") else 0,
                    1 if site["coordination_number"] == 0 else 0, ligands))
        return entry_id, rows, None
    except Exception as exc:  # per-file failure must not kill the batch
        entry_id = os.path.basename(path).replace(".cif.gz", "")
        return entry_id, [], f"{type(exc).__name__}: {exc}"[:150]


def build_index(list_json: str = LIST_JSON, download_dir: str = DOWNLOAD_DIR,
                db_path: str = DB_PATH, workers: int = 4,
                limit: int | None = None, commit_every: int = 25) -> dict:
    """Stage 3: multiprocess QC over the downloaded subset -> SQLite index.

    Resumable: entries already in the DB are skipped. On Windows this uses
    spawn - call it via the CLI (main) rather than from a notebook."""
    mapping = _load_mapping(list_json) if os.path.exists(list_json) else {}
    metals = sorted({m for m in mapping} ) or list(DEFAULT_METALS)
    files = sorted(f for f in os.listdir(download_dir)
                   if f.endswith((".cif.gz", ".cif")))
    paths = [os.path.join(download_dir, f) for f in files]
    if limit:
        paths = paths[:limit]
    con = _init_db(db_path)
    done = {r[0] for r in con.execute("SELECT DISTINCT entry_id FROM sites")}
    todo = [p for p in paths
            if os.path.basename(p).replace(".cif.gz", "").replace(".cif", "")
            not in done]
    print(f"build: {len(todo)} of {len(paths)} files to index "
          f"({len(done)} already done), metals={metals}, {workers} workers")
    n_rows = 0
    errors: list[dict] = []
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_process_file, p, metals) for p in todo]
        for i, fut in enumerate(futures):
            entry_id, rows, err = fut.result()
            if err:
                errors.append({"entry": entry_id, "error": err})
            elif rows:
                con.executemany(INSERT_SQL, rows)
                n_rows += len(rows)
            if (i + 1) % commit_every == 0:
                con.commit()
                print(f"  {i + 1}/{len(todo)} files, {n_rows} sites so far")
    con.commit()
    total = con.execute("SELECT COUNT(*) FROM sites").fetchone()[0]
    con.close()
    out = {"files_scanned": len(todo), "new_rows": n_rows,
           "total_rows": total, "errors": errors[:10]}
    print(f"build done: {out['new_rows']} new sites, {total} total")
    return out
def status(list_json: str = LIST_JSON, download_dir: str = DOWNLOAD_DIR,
           db_path: str = DB_PATH) -> dict[str, Any]:
    """Inventory: what the list/sync/build stages have produced so far."""
    out: dict[str, Any] = {}
    if os.path.exists(list_json):
        mapping = _load_mapping(list_json)
        out["list"] = {m: len(v) for m, v in mapping.items()}
        out["unique_entries"] = len({i for v in mapping.values() for i in v})
    files = [f for f in os.listdir(download_dir)
             if f.endswith((".cif.gz", ".cif"))] if os.path.isdir(download_dir) else []
    out["downloaded"] = {"n_files": len(files),
                         "total_mb": round(sum(
                             os.path.getsize(os.path.join(download_dir, f))
                             for f in files) / 1e6, 1)}
    if os.path.exists(db_path):
        con = sqlite3.connect(db_path)
        out["db_rows"] = con.execute("SELECT COUNT(*) FROM sites").fetchone()[0]
        out["db_entries"] = con.execute(
            "SELECT COUNT(DISTINCT entry_id) FROM sites").fetchone()[0]
        out["per_metal"] = dict(con.execute(
            "SELECT metal, COUNT(*) FROM sites GROUP BY metal").fetchall())
        out["geometry_distribution"] = dict(con.execute(
            "SELECT geometry, COUNT(*) FROM sites GROUP BY geometry "
            "ORDER BY COUNT(*) DESC").fetchall())
        con.close()
        out["db_path"] = db_path
    return out


def query(metal: str, geometry: str | None = None, min_cn: int = 2,
          max_cn: int = 12, limit: int = 20,
          include_artifacts: bool = False,
          db_path: str = DB_PATH) -> dict[str, Any]:
    """Query the index: e.g. all octahedral ZN sites with 4 <= CN <= 6.

    Flag-don't-fail defaults: CN-0 rows are classified artifact=1
    (isolated-ion artifact) and excluded unless include_artifacts=True;
    min_cn defaults to 2 so CN-0/CN-1 rows are hidden by default while
    CN-1 unusual chemistry stays queryable via min_cn=1.
    """
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"index not built yet: {db_path}")
    con = sqlite3.connect(db_path)
    q = ("SELECT entry_id, metal, site, coordination_number, geometry, "
         "rmsd_deg, n_imidazole, n_carboxylate, n_water, n_thiolate, "
         "n_adjacent_metal, bidentate, artifact, ligands "
         "FROM sites WHERE metal = ?")
    args: list = [metal.upper()]
    if include_artifacts:
        q += " AND (coordination_number BETWEEN ? AND ? OR artifact = 1)"
    else:
        q += " AND coordination_number BETWEEN ? AND ? AND artifact = 0"
    args += [min_cn, max_cn]
    if geometry:
        q += " AND geometry = ?"
        args.append(geometry.lower())
    q += " LIMIT ?"
    args.append(limit)
    rows = con.execute(q, args).fetchall()
    con.close()
    cols = ["entry_id", "metal", "site", "cn", "geometry", "rmsd_deg",
            "n_imidazole", "n_carboxylate", "n_water", "n_thiolate",
            "n_adjacent_metal", "bidentate", "artifact", "ligands"]
    return {"n_results": len(rows), "results": [dict(zip(cols, r))
                                                for r in rows]}


def main() -> None:
    ap = argparse.ArgumentParser(
        description="PDB-wide metal-site index builder "
                    "(list / sync / build / status / query)")
    ap.add_argument("stage", choices=["list", "sync", "build", "status",
                                      "query"])
    ap.add_argument("--metals", default=",".join(DEFAULT_METALS),
                    help="comma-separated metal symbols (list stage)")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap per metal (list) or total entries (sync/build)")
    ap.add_argument("--workers", type=int, default=None,
                    help="parallel workers (sync: threads, build: processes)")
    ap.add_argument("--max-mb", type=float, default=50.0,
                    help="skip downloads larger than this (sync)")
    ap.add_argument("--metal", default="ZN", help="query stage: metal symbol")
    ap.add_argument("--geometry", default=None,
                    help="query stage: geometry name (e.g. octahedral)")
    ap.add_argument("--min-cn", type=int, default=2)
    ap.add_argument("--max-cn", type=int, default=12)
    ap.add_argument("--include-artifacts", action="store_true",
                    help="query stage: also return CN-0 artifact rows")
    args = ap.parse_args()

    if args.stage == "list":
        fetch_metal_entry_list([m.strip().upper() for m in
                                args.metals.split(",") if m.strip()],
                               limit_per_metal=args.limit)
    elif args.stage == "sync":
        sync_entries(workers=args.workers or 8, max_mb=args.max_mb,
                     limit=args.limit)
    elif args.stage == "build":
        build_index(workers=args.workers or 4, limit=args.limit)
    elif args.stage == "status":
        print(json.dumps(status(), indent=2))
    elif args.stage == "query":
        print(json.dumps(query(args.metal, geometry=args.geometry,
                               min_cn=args.min_cn, max_cn=args.max_cn,
                               include_artifacts=args.include_artifacts,
                               limit=args.limit or 20), indent=2))


if __name__ == "__main__":
    main()



