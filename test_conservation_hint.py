import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import premium_tools

# 1. Canned UniProt JSON: all six wanted feature types
raw = {
    "primaryAccession": "P00533",
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Epidermal growth factor receptor"}}},
    "organism": {"scientificName": "Homo sapiens"},
    "sequence": {"length": 1210, "value": "M" * 1210},
    "comments": [{"texts": ["Receptor tyrosine kinase."]}],
    "features": [
        {"type": "Metal binding", "description": "Magnesium",
         "location": {"start": {"value": 100}, "end": {"value": 100}}},
        {"type": "Active site", "description": "Proton acceptor",
         "location": {"start": {"value": 200}, "end": {"value": 200}}},
        {"type": "Binding site", "description": "ATP",
         "location": {"start": {"value": 300}, "end": {"value": 340}}},
        {"type": "Disulfide bond", "description": "By similarity",
         "location": {"start": {"value": 400}, "end": {"value": 450}}},
        {"type": "Mutagenesis", "description": "Loss of activity",
         "location": {"start": {"value": 500}, "end": {"value": 500}}},
        {"type": "Natural variant", "description": "In cancer",
         "location": {"start": {"value": 600}, "end": {"value": 600}}},
        {"type": "Signal peptide", "description": "ignored type",
         "location": {"start": {"value": 1}, "end": {"value": 20}}},
    ],
}
s = premium_tools.uniprot_summary(raw)
assert s["accession"] == "P00533", s["accession"]
print("1. conservation_hint summarizer accession: OK")
assert len(s["metal_binding_sites"]) == 1, s["metal_binding_sites"]
print("2. conservation_hint legacy metal_binding_sites intact: OK")
assert set(s["features_by_type"]) == {"Metal binding", "Active site",
    "Binding site", "Disulfide bond", "Mutagenesis", "Natural variant"}, \
    set(s["features_by_type"])
print("3. conservation_hint features_by_type (6 types, Signal ignored): OK")
assert s["features_by_type"]["Binding site"][0]["end"] == 340
print("4. conservation_hint start/end positions: OK")

# 5. Premium gate: conservation_hint without a key -> structured error
import os
os.environ.pop("PROTEIN_TOOLS_API_KEY", None)
os.environ.pop("PROTEIN_TOOLS_RELAY_URL", None)
# import the tool function without starting the MCP server: load module
# functions directly (server imports mcp SDK — install it via requirements).
import server as srv
import json
out = json.loads(srv.conservation_hint("P00533"))
assert "error" in out and "premium" in out["error"].lower(), out
print("5. conservation_hint premium gate (no key -> error): OK")

print()
print("ALL conservation_hint TESTS PASSED")