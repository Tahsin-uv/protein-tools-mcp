"""protein-tools relay — hosted premium backend (Phase A MVP, manual keys).

Deploys on Modal as a FastAPI service. This is the thing that gets sold:
- Hosts the PDB-wide metal-site index (metal-index-vol, 38,143 entries)
- Validates API keys (SHA-256 hashed; stored on protein-tools-store volume)
- Enforces monthly quotas per tier + a simple per-minute rate limit
- Meters every call (usage visible via GET /v1/usage)
- Server-side 24 h cache for UniProt / AlphaFold / RCSB upstreams

Usage:
  modal deploy relay.py                                        # deploy API
  modal run relay.py --stage issue --name beta1 --tier beta    # issue key
  modal run relay.py --stage list                              # inventory
  modal run relay.py --stage revoke --key pt_live_xxx          # revoke

NOTE: do NOT add `from __future__ import annotations` to this file —
FastAPI must resolve `request: Request` at def-time (closure-local import).
"""

import hashlib
import json
import os
import secrets
import sqlite3
import time
from typing import Any

import modal

KEYS_VOL = modal.Volume.from_name("protein-tools-store", create_if_missing=True)
IDX_VOL = modal.Volume.from_name("metal-index-vol")

STORE = "/store"
CACHE = f"{STORE}/cache"
DB = f"{STORE}/keys.db"
INDEX_DB = "/vol/merged/metal_sites.db"

image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "fastapi").add_local_python_source(
    "builder_metal_index", "cif_tools", "pdb_utils", "premium_tools")

app = modal.App("protein-tools-relay", image=image)

# ----------------------------- pure logic -----------------------------------

TIERS = {"beta": 1000, "pro": 10000, "lab": 50000}
RATE_LIMIT_PER_MIN = 60
UPSTREAM_TTL_H = 24.0
UA = "protein-tools-relay/0.1"


def _hash_key(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def _gen_key() -> str:
    return "pt_live_" + secrets.token_hex(16)


def _period() -> str:
    return time.strftime("%Y-%m", time.gmtime())


def _quota_for(tier: str) -> int:
    return TIERS.get(tier, TIERS["beta"])


def _effective_usage(row: sqlite3.Row) -> dict[str, Any]:
    """Quota state for a key row, resetting the monthly counter when the
    period rolled over."""
    period = _period()
    used, row_period = row["calls_used"], row["period"]
    if row_period != period:
        used = 0
        con = _connect()
        con.execute("UPDATE keys SET calls_used = 0, period = ? "
                    "WHERE key_hash = ?", (period, row["key_hash"]))
        con.commit()
        con.close()
    return {"calls_used": used, "quota": _quota_for(row["tier"]),
            "period": period}


def _connect() -> sqlite3.Connection:
    con = sqlite3.connect(DB, timeout=30)
    con.row_factory = sqlite3.Row
    con.execute("CREATE TABLE IF NOT EXISTS keys ("
                "key_hash TEXT PRIMARY KEY, name TEXT, tier TEXT, "
                "created REAL, revoked INTEGER DEFAULT 0, "
                "calls_used INTEGER DEFAULT 0, period TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS calls ("
                "id INTEGER PRIMARY KEY AUTOINCREMENT, key_hash TEXT, "
                "ts REAL, endpoint TEXT, status INTEGER)")
    return con


def _init_db() -> None:
    os.makedirs(os.path.dirname(DB), exist_ok=True)
    con = _connect()
    con.commit()
    con.close()


def issue_key(name: str, tier: str = "beta") -> dict[str, Any]:
    """Generate a new key, store only its hash. Plaintext shown once."""
    _init_db()
    key = _gen_key()
    con = _connect()
    con.execute("INSERT INTO keys (key_hash, name, tier, created, revoked, "
                "calls_used, period) VALUES (?, ?, ?, ?, 0, 0, ?)",
                (_hash_key(key), name, tier, time.time(), _period()))
    con.commit()
    con.close()
    return {"key": key, "name": name, "tier": tier,
            "quota_per_month": _quota_for(tier)}


def revoke_key(key: str) -> dict[str, Any]:
    h = _hash_key(key)
    con = _connect()
    cur = con.execute("UPDATE keys SET revoked = 1 WHERE key_hash = ?", (h,))
    con.commit()
    n = cur.rowcount
    con.close()
    return {"revoked": bool(n)}


def list_keys() -> list[dict[str, Any]]:
    _init_db()
    con = _connect()
    rows = con.execute("SELECT key_hash, name, tier, created, revoked, "
                       "calls_used, period FROM keys ORDER BY created "
                       "DESC").fetchall()
    con.close()
    return [{"key_hash": r["key_hash"][:12] + "...", "name": r["name"],
             "tier": r["tier"], "revoked": bool(r["revoked"]),
             "calls_used": r["calls_used"], "period": r["period"]}
            for r in rows]


def get_valid_key(key: str) -> dict[str, Any]:
    """Validate a bearer key; raises PermissionError when invalid, revoked,
    or over quota."""
    _init_db()
    h = _hash_key(key)
    con = _connect()
    row = con.execute("SELECT * FROM keys WHERE key_hash = ?", (h,)).fetchone()
    con.close()
    if row is None:
        raise PermissionError("unknown API key")
    if row["revoked"]:
        raise PermissionError("this API key has been revoked")
    usage = _effective_usage(row)
    if usage["calls_used"] >= usage["quota"]:
        raise PermissionError(
            f"monthly quota exhausted ({usage['quota']}/month on the "
            f"'{row['tier']}' tier) - upgrade or wait for the next period")
    return {"key_hash": h, "name": row["name"], "tier": row["tier"],
            **usage}


# ---------------------------- upstream cache --------------------------------

def _cached_get(url: str) -> bytes:
    """Server-side disk cache (shared across ALL users = quota-efficient)."""
    import urllib.request
    os.makedirs(CACHE, exist_ok=True)
    key = hashlib.md5(url.encode("utf-8")).hexdigest()
    blob = os.path.join(CACHE, key + ".cache")
    meta = blob + ".meta"
    if os.path.exists(blob) and os.path.exists(meta):
        try:
            if time.time() - float(open(meta).read().strip()) \
                    < UPSTREAM_TTL_H * 3600:
                return open(blob, "rb").read()
        except (OSError, ValueError):
            pass
    last: Exception | None = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            with open(blob, "wb") as f:
                f.write(data)
            with open(meta, "w") as f:
                f.write(str(time.time()))
            return data
        except Exception as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise last


def meter(key_hash: str, endpoint: str, status: int) -> None:
    con = _connect()
    con.execute("INSERT INTO calls (key_hash, ts, endpoint, status) "
                "VALUES (?, ?, ?, ?)", (key_hash, time.time(), endpoint,
                                        status))
    con.execute("UPDATE keys SET calls_used = calls_used + 1 "
                "WHERE key_hash = ?", (key_hash,))
    con.commit()
    con.close()


# ------------------------------- web app ------------------------------------

@app.function(image=image, volumes={STORE: KEYS_VOL, "/vol": IDX_VOL},
               timeout=600)
@modal.asgi_app(label="relay")
def api() -> Any:
    from fastapi import FastAPI, HTTPException, Request
    from fastapi.responses import JSONResponse
    import builder_metal_index as metal_index
    import pdb_utils
    import premium_tools as pt

    web = FastAPI(title="protein-tools relay", version="0.1.0")
    rate_window: dict[str, list[float]] = {}

    def auth(request: Request) -> dict[str, Any]:
        authz = request.headers.get("authorization", "")
        key = authz.removeprefix("Bearer ").strip()
        if not key:
            raise HTTPException(401, "missing Authorization: Bearer <key>")
        now = time.time()
        win = [t for t in rate_window.get(key[:12], []) if now - t < 60.0]
        if len(win) >= RATE_LIMIT_PER_MIN:
            raise HTTPException(429, "rate limit: 60 requests/minute")
        win.append(now)
        rate_window[key[:12]] = win
        try:
            info = get_valid_key(key)
        except PermissionError as exc:
            raise HTTPException(403, str(exc)) from exc
        return info

    def guarded(request: Request, endpoint: str, fn) -> JSONResponse:
        info = auth(request)
        try:
            payload, status = fn()
        except HTTPException:
            raise
        except Exception as exc:
            payload = {"error": f"{type(exc).__name__}: {exc}",
                       "hint": "retry later or report this message"}
            status = 502
        try:
            meter(info["key_hash"], endpoint, status)
        except Exception:
            pass  # metering must never break a served response
        return JSONResponse(payload, status_code=status)

    @web.get("/v1/health")
    def health():
        return {"ok": True, "service": "protein-tools-relay",
                "version": "0.1.0"}

    @web.get("/v1/usage")
    def usage(request: Request):
        info = auth(request)
        meter(info["key_hash"], "usage", 200)
        return {"name": info["name"], "tier": info["tier"],
                "calls_used": info["calls_used"], "quota": info["quota"],
                "period": info["period"]}

    @web.get("/v1/query_metal_sites")
    def query_metal_sites(request: Request, metal: str = "ZN",
                          geometry: str = "", min_cn: int = 2,
                          max_cn: int = 12, limit: int = 20,
                          include_artifacts: bool = False):
        def run():
            r = metal_index.query(
                metal, geometry=geometry or None, min_cn=min_cn,
                max_cn=max_cn, limit=min(limit, 200),
                include_artifacts=include_artifacts, db_path=INDEX_DB)
            return r, 200
        return guarded(request, "query_metal_sites", run)

    @web.get("/v1/index_stats")
    def index_stats(request: Request):
        def run():
            return metal_index.status(db_path=INDEX_DB), 200
        return guarded(request, "index_stats", run)

    @web.get("/v1/uniprot/{accession}")
    def uniprot(request: Request, accession: str):
        def run():
            d = json.loads(_cached_get(
                f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
            return pt.uniprot_summary(d), 200
        return guarded(request, "uniprot", run)

    @web.get("/v1/alphafold/{accession}")
    def alphafold(request: Request, accession: str):
        def run():
            text = _cached_get(
                f"https://alphafold.ebi.ac.uk/files/AF-{accession}"
                f"-F1-model_v4.pdb").decode(errors="replace")
            return pdb_utils.plddt_summary(text), 200
        return guarded(request, "alphafold", run)

    @web.get("/v1/rcsb/{pdb_id}")
    def rcsb(request: Request, pdb_id: str):
        def run():
            d = json.loads(_cached_get(
                f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"))
            return pt.rcsb_summary(d), 200
        return guarded(request, "rcsb", run)

    return web


# ---------------------------- CLI management --------------------------------

@app.function(image=image, volumes={STORE: KEYS_VOL}, timeout=300)
def issue_key_remote(name: str, tier: str) -> dict[str, Any]:
    return issue_key(name, tier)


@app.function(image=image, volumes={STORE: KEYS_VOL}, timeout=300)
def revoke_key_remote(key: str) -> dict[str, Any]:
    return revoke_key(key)


@app.function(image=image, volumes={STORE: KEYS_VOL}, timeout=300)
def list_keys_remote() -> list[dict[str, Any]]:
    return list_keys()


@app.local_entrypoint()
def main(stage: str = "help", name: str = "", tier: str = "beta",
         key: str = "") -> None:
    if stage == "issue":
        print(json.dumps(issue_key_remote.remote(name or "unnamed", tier),
                         indent=1))
    elif stage == "list":
        print(json.dumps(list_keys_remote.remote(), indent=1))
    elif stage == "revoke":
        if not key:
            raise SystemExit("--key is required for revoke")
        print(json.dumps(revoke_key_remote.remote(key), indent=1))
    else:
        print("stages: issue (--name, --tier beta|pro|lab) | list | revoke "
              "(--key) - then: modal deploy relay.py")
