"""Protein Lab — plain-English analysis for students (no MCP setup needed).

Wraps the tested analysis modules (seq_utils, pdb_utils, cif_tools) in a
Streamlit UI: paste/upload -> traffic-light verdict cards -> "What does this
mean?" explainers. Demo datasets are the synthetic fixtures/ files, so
classroom mode works fully offline.

Local run:   streamlit run protein_lab.py
HF Spaces:   this file IS the Spaces entry point (see README section).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import lab_seq  # noqa: E402
import lab_structure  # noqa: E402

try:
    import streamlit as st
except ImportError:  # pragma: no cover - CLI fallback hint
    print("Protein Lab needs streamlit: pip install streamlit")
    raise SystemExit(1)

st.set_page_config(page_title="Protein Lab", page_icon="🧬",
                   layout="centered")

st.title("🧬 Protein Lab")
st.caption("Upload a file (or try the demo data) — get a plain-English "
           "answer. No installs, no code, no JSON.")

tab_seq, tab_struct = st.tabs(["🔤 Sequence", "🏗️ Structure"])

with tab_seq:
    lab_seq.render()

with tab_struct:
    lab_structure.render()

st.divider()
st.caption("Protein Lab — the web UI for "
           "[protein-tools](https://github.com/tahsinhassan18/protein-tools-mcp). "
           "Analysis by the free core; figures use your data only — "
           "nothing is uploaded anywhere.")
