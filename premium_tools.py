"""Shared summarizers for premium lookups (used by both server.py and the
Modal relay, so client-side and relay-side responses never drift)."""
from __future__ import annotations

from typing import Any


def uniprot_summary(d: dict[str, Any]) -> dict[str, Any]:
    """Reduce a raw UniProtKB JSON entry to the fields the MCP tool exposes."""
    desc = d.get("proteinDescription", {})
    name = (desc.get("recommendedName") or desc.get("submissionNames") or [{}])
    name = name.get("fullName", {}).get("value") if isinstance(name, dict) \
        else name[0].get("fullName", {}).get("value")
    metals = [{"description": f.get("description"),
               "position": f.get("location", {}).get("start", {}).get("value")}
              for f in d.get("features", []) if f.get("type") == "Metal binding"]
    # Grouped feature table for conservation_hint (additive: existing
    # callers using metal_binding_sites are unaffected).
    wanted = {"Metal binding", "Active site", "Binding site",
              "Disulfide bond", "Mutagenesis", "Natural variant"}
    by_type: dict[str, list] = {}
    for f in d.get("features", []):
        t = f.get("type")
        if t not in wanted:
            continue
        loc = f.get("location", {})
        by_type.setdefault(t, []).append({
            "description": f.get("description"),
            "start": (loc.get("start") or {}).get("value"),
            "end": (loc.get("end") or {}).get("value"),
        })
    seq = d.get("sequence", {})
    return {
        "accession": d.get("primaryAccession"),
        "name": name,
        "organism": d.get("organism", {}).get("scientificName"),
        "length": seq.get("length"),
        "sequence_head": (seq.get("value") or "")[:80] + "...",
        "function": (d.get("comments") or [{}])[0].get("texts", [""])[0][:400],
        "metal_binding_sites": metals,
        "features_by_type": by_type,
    }


def rcsb_summary(d: dict[str, Any]) -> dict[str, Any]:
    """Reduce a raw RCSB entry JSON to the fields the MCP tool exposes."""
    info = d.get("rcsb_entry_info", {})
    res = info.get("resolution_combined") or [None]
    return {
        "pdb_id": d.get("rcsb_id"),
        "title": d.get("struct", {}).get("title"),
        "method": (d.get("exptl") or [{}])[0].get("method"),
        "resolution_A": res[0],
        "deposited": d.get("rcsb_accession_info", {}).get("deposit_date"),
    }
