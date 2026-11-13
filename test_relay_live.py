"""Live verification of the deployed relay (Phase A acceptance test).
Run only after: modal deploy relay.py + key issue."""
import json
import os
import sys
import urllib.error
import urllib.request

BASE = "https://tahsinhassanakash--relay.modal.run"
KEY = "pt_live_4e051bc4513e28647b504bbb983b2101"


def get(path: str, key: str | None = KEY):
    req = urllib.request.Request(BASE + path)
    if key:
        req.add_header("Authorization", f"Bearer {key}")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read().decode())
    except urllib.error.HTTPError as exc:
        try:
            return exc.code, json.loads(exc.read().decode())
        except Exception:
            return exc.code, {}


# 1. health: public, no auth
s, d = get("/v1/health", key=None)
assert s == 200 and d["ok"] is True, (s, d)
print("1 health: OK", d)

# 2. no key -> 401
s, d = get("/v1/query_metal_sites?metal=ZN", key=None)
assert s == 401, (s, d)
print("2 no-key -> 401: OK")

# 3. bad key -> 403
s, d = get("/v1/query_metal_sites?metal=ZN", key="pt_live_bogus")
assert s == 403, (s, d)
print("3 bad-key -> 403: OK", d.get("error", ""))

# 4. good key -> flagship query against the hosted index
s, d = get("/v1/query_metal_sites?metal=ZN&geometry=tetrahedral&min_cn=4&max_cn=4&limit=3")
assert s == 200 and d["n_results"] == 3, (s, d)
assert d["results"][0]["entry_id"] == "10AJ"
assert d["results"][0]["n_thiolate"] == 4
print("4 index query: OK", json.dumps(d["results"][0])[:120])

# 5. index stats
s, d = get("/v1/index_stats")
assert s == 200 and d.get("db_rows") == 1097659, (s, d)
print("5 index stats: OK rows =", d["db_rows"])

# 6. usage metering visible
s, d = get("/v1/usage")
assert s == 200 and d["calls_used"] >= 2, (s, d)
print("6 usage: OK", d)

# 7. premium upstream via relay (cached server-side)
s, d = get("/v1/uniprot/A8ABE7")
assert s == 200 and d["accession"] == "A8ABE7", (s, d)
print("7 uniprot via relay: OK", d["organism"])

print("RELAY-LIVE OK: all acceptance checks passed")
