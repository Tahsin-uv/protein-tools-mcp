"""Unit tests for the relay's pure logic + shared premium summarizers.
Runs fully offline (key DB redirected to a temp dir)."""
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import relay  # noqa: E402
import premium_tools  # noqa: E402

tmp = tempfile.mkdtemp()
relay.DB = os.path.join(tmp, "keys.db")
relay.CACHE = os.path.join(tmp, "cache")

# --- key lifecycle -----------------------------------------------------------
assert relay._hash_key("abc") == relay._hash_key("abc")
assert relay._hash_key("abc") != relay._hash_key("abd")

k = relay.issue_key("unit-t1", "beta")
key = k["key"]
assert key.startswith("pt_live_") and k["quota_per_month"] == 1000

info = relay.get_valid_key(key)
assert info["tier"] == "beta" and info["calls_used"] == 0

relay.meter(info["key_hash"], "query_metal_sites", 200)
info2 = relay.get_valid_key(key)
assert info2["calls_used"] == 1, info2

try:
    relay.get_valid_key("pt_live_unknown")
    raise AssertionError("unknown key should fail")
except PermissionError:
    pass

relay.revoke_key(key)
try:
    relay.get_valid_key(key)
    raise AssertionError("revoked key should fail")
except PermissionError:
    pass

# --- tiers -------------------------------------------------------------------
assert relay._quota_for("pro") == 10000
assert relay._quota_for("lab") == 50000
assert relay._quota_for("nonsense") == 1000  # safe default

# --- monthly period reset ----------------------------------------------------
k2 = relay.issue_key("unit-t2", "pro")["key"]
h2 = relay._hash_key(k2)
con = relay._connect()
con.execute("UPDATE keys SET calls_used = 99999, period = '1999-01' "
            "WHERE key_hash = ?", (h2,))
con.commit()
con.close()
info3 = relay.get_valid_key(k2)
assert info3["calls_used"] == 0, info3  # old period -> counter reset

# --- shared summarizers ------------------------------------------------------
u = premium_tools.uniprot_summary({
    "primaryAccession": "P12345",
    "organism": {"scientificName": "Test org"},
    "proteinDescription": {"recommendedName": {"fullName": {"value": "X"}}},
    "sequence": {"length": 10, "value": "ABCDEFGHIJ"},
    "features": [{"type": "Metal binding", "description": "ZN",
                  "location": {"start": {"value": 7}}}],
})
assert u["accession"] == "P12345" and u["length"] == 10
assert u["metal_binding_sites"] == [{"description": "ZN", "position": 7}]

r = premium_tools.rcsb_summary({
    "rcsb_id": "1ABC", "struct": {"title": "T"},
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "rcsb_entry_info": {"resolution_combined": [1.5]},
    "rcsb_accession_info": {"deposit_date": "2020-01-01"},
})
assert r == {"pdb_id": "1ABC", "title": "T",
             "method": "X-RAY DIFFRACTION", "resolution_A": 1.5,
             "deposited": "2020-01-01"}, r

print("RELAY-UNIT OK: keys, quota reset, tiers, summarizers")
