"""Sequence tab for Protein Lab: paste/type/upload -> plain-English cards."""
from __future__ import annotations

import os

import seq_utils

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")

_TYPE_WORDS = {
    "protein": ("a **protein** (a chain of amino acids)", "🟢"),
    "dna": ("**DNA** (the genetic code molecule)", "🟢"),
    "rna": ("**RNA** (genetic messenger molecule)", "🟢"),
    "ambiguous": ("**unclear** — it could be DNA or RNA", "🟡"),
    "unknown": ("**not readable** — check the input", "🔴"),
}

_GLOSSARY = {
    "Amino acid": "One building block of a protein. 20 types exist; "
                  "their order (the sequence) decides what the protein does.",
    "Charge at pH 7": "Whether the molecule is overall positive or "
                      "negative in a cell-like solution.",
    "GRAVY": "Positive means water-hating (oily, sits in membranes); "
             "negative means water-loving (dissolves in cell fluid).",
    "pI": "The acidity (pH) where the protein has zero net charge — "
          "useful for purification in the lab.",
    "Molecular weight": "How heavy one copy of the molecule is, in "
                        "Daltons (~110 Da per amino acid on average).",
}


def _demo_text() -> str:
    p = os.path.join(FIX, "sample.fasta")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return f.read()
    return ">demo_protein\nMKLIVFWYHEPQDNRSACDE\n"

def render() -> None:
    import streamlit as st

    # Theme selector (Streamlit built-in)
    themes = ["light", "dark", "auto"]
    theme = st.sidebar.radio("App theme", themes,
                          index=["light", "dark", "auto"].index(st.get_option("theme.base")))
    st.set_page_config(page_title="Protein Lab", page_icon="🧬", layout="wide",
                      initial_sidebar_state="expanded")

    # Session history: keep last 3 inputs
    if "lab_seq_history" not in st.session_state:
        st.session_state.lab_seq_history = []

    st.sidebar.header("History")
    if st.session_state.lab_seq_history:
        for i, hist in enumerate(reversed(st.session_state.lab_seq_history[-3:])):
            if st.sidebar.button(f"Recall #{len(st.session_state.lab_seq_history)-i}", key=f"recall_{i}"):
                st.session_state._lab_seq = hist
                st.rerun()
    if st.sidebar.button("Clear history"):
        st.session_state.lab_seq_history = []
        st.rerun()

    st.header("What is my sequence?")
    st.write("Paste letters (like `MKLIVFWYHEPQDNRS`), a FASTA record "
             "(starting with `>`), or upload a `.fasta` file.")

    if "_lab_seq" not in st.session_state:
        st.session_state["_lab_seq"] = ""

    # Handle actions that modify session state BEFORE creating the widget,
    # so the widget picks up the new value on rerun.
    demo_clicked = st.button("Try demo data instead")
    if demo_clicked:
        st.session_state["_lab_seq"] = _demo_text()
        st.session_state.lab_seq_history.append(st.session_state._lab_seq)
        # Keep only last 10
        st.session_state.lab_seq_history = st.session_state.lab_seq_history[-10:]

    up = st.file_uploader("...or upload a .fasta file",
                          type=["fasta", "fa", "txt"])
    if up is not None:
        st.session_state["_lab_seq"] = up.read().decode("utf-8", errors="replace")
        st.session_state.lab_seq_history.append(st.session_state._lab_seq)
        st.session_state.lab_seq_history = st.session_state.lab_seq_history[-10:]

    src = st.text_area("Sequence or FASTA", height=120,
                       placeholder=">my_protein\nMKLIVFWYHEPQDNRSACDE",
                       key="_lab_seq")

    submitted = st.button("Analyze Sequence")

    src = st.session_state.get("_lab_seq", "")
    if not src.strip() and not submitted:
        st.info("Paste a sequence, upload a file, or try the demo data.")
        return

    text = src.strip()
    try:
        if text.startswith(">") or "\n>" in text:
            recs = seq_utils.parse_fasta(text, is_file=False)
            st.subheader(f"{len(recs)} sequence(s) found")
            for r in recs:
                _render_record(r["id"], r["sequence"],
                               r["type"]["type"])
        else:
            seq = seq_utils.normalize_sequence(text)
            label, detail = seq_utils.classify_sequence_type(seq)
            if label == "protein":
                _render_record("your input", seq, label)
            else:
                words, dot = _TYPE_WORDS.get(label, (label, ""))
                st.markdown(f"### {dot} This looks like {words}.")
                st.caption(detail.note)
    except (ValueError, FileNotFoundError) as exc:
        st.error(f"Couldn't read that input: {exc}")
        st.caption("Check for stray spaces, numbers, or a missing `>` "
                   "header line in FASTA input.")

    # Format guide: what to paste/upload
    with st.expander("File format guide"):
        st.markdown("""
        **Accepted formats**:
        - **FASTA**: Start with `>` followed by a header line, then the sequence on the next line(s).
          Example: `>sp|P68871|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=1`
        - **Plain text**: Just the letter code, e.g. `MKLIVFWYHEPQDNRSACDE`
        - **PDB/mmCIF**: Upload `.pdb` or `.cif` files in the Structure tab.
          Must have `ATOM` or `REMARK` records, and residue names must be standard 3-letter codes.

        **Tips**:
        - Remove any leading/trailing spaces or numbers before the header.
        - FASTA headers with `|` or spaces are fine — only the first `>` matters.
        - If your sequence looks like "amino a..." you likely have a truncation issue —
          the app now shows "aa" instead of truncating.
        - For PDB files: if you see "Format errors", check that the file has proper
          `HEADER`, `TITLE`, `COMPND`, and `SOURCE` records for best results.
        """)

    with st.expander("Glossary - what do these words mean?"):
        for term, meaning in _GLOSSARY.items():
            st.markdown(f"**{term}:** {meaning}")

# New: Detailed help/tooltips for each metric
with st.expander("Help & interpretation"):
    st.markdown(f"""
    **Length**: Number of amino acids in your protein.  
    **Weight**: Molecular weight in Daltons.  ~110 Da per residue is the average.  
    **Charge at pH 7**: Net charge the protein would have in a typical cellular environment (pH ~7.4).  
    **GRAVY**: Grand Average of Hydropathy. Positive = water-hating (oily), negative = water-loving.  
      Used to predict transmembrane regions — if strongly positive, the protein may embed in membranes.  
    **pI (isoelectric point)**: The pH at which the protein has zero net charge.  
      Below the pI the protein is positively charged; above it is negatively charged.  
      Important for chromatographic purification — proteins precipitate near their pI.  
    **ε₂₈₀**: Molar extinction coefficient at 280 nm. Determined by Trp (×5500), Tyr (×1490), and Cys (×125 if oxidized).  
      Used to estimate protein concentration from A₂₈₀ measurements.  
    **A₂₈₀ @ 0.1%**: Absorbance of a 1 mg/mL solution in a 1 cm cell.  
      Rough check: if ε₂₈₀ is correct, A₂₈₀ ≈ 1.0 for 1 mg/mL proteins with typical Trp/Tyr content.  
    **Amino acid types**: Hydrophobic (A,V,L,I,M,F,W), Polar (N,C,Q,Y,S,T), Acidic (D,E), Basic (K,R,H),  
      Aromatic (F,W,Y), Special (C, glycine, proline). The chart shows your distribution.  
    **Class fractions**: The dominant class tells you what your protein "likes" — e.g. hydrophobic proteins  
      tend to stick to membranes, acidic proteins bind cations, basic proteins bind nucleic acids.  
    **Charge direction**: If the protein is positively charged, it will tend to associate with  
      negatively charged surfaces (like DNA or phospholipid heads); if negative, the opposite.
    """)


def _render_record(rec_id: str, seq: str, label: str) -> None:
    import streamlit as st

    words, dot = _TYPE_WORDS.get(label, (label, ""))
    st.markdown(f"### {dot} `{rec_id}` looks like {words}.")
    st.caption(f"{len(seq)} letters long.")
    if label != "protein":
        return
    try:
        c = seq_utils.aa_composition(seq)
    except ValueError as exc:
        st.warning(str(exc))
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Length", f"{c['length']} aa")
    c2.metric("Weight", f"{c['molecular_weight_da']:,.0f} Da")
    ch = c.get("charge_at_pH7_estimate")
    c3.metric("Charge at pH 7", f"{ch:+.1f}" if ch is not None else "-")
    g1, g2 = st.columns(2)
    gravy = c.get("gravy")
    g1.metric("GRAVY (water-hating <-> water-loving)",
              f"{gravy:+.2f}" if gravy is not None else "-",
              help=_GLOSSARY["GRAVY"])
    pi = c.get("pi_estimate")
    g2.metric("pI (zero-charge acidity)", f"{pi:.2f}" if pi else "-",
              help=_GLOSSARY["pI"])

    # Extinction coefficient (Pace et al., 1995)
    try:
        ec = seq_utils.extinction_coefficient(seq)
        e1, e2 = st.columns(2)
        e1.metric("ε₂₈₀ (M⁻¹ cm⁻¹)", f"{ec['extinction_coefficient_m1cm1']:,.0f}",
                  help="Molar extinction coefficient at 280 nm")
        e2.metric("A₂₈₀ @ 0.1% (1 mg/mL)", f"{ec['absorbance_0_1_percent']:.2f}",
                  help="Absorbance of 1 mg/mL solution in 1 cm cell")
        st.caption(f"  Trp: {ec['n_trp']} | Tyr: {ec['n_tyr']} | "
                   f"Cys (reduced): {ec['n_cys']} — {ec['note']}")
    except ValueError as exc:
        st.caption(f"Extinction coefficient: skipped ({exc})")
    st.write("**Building blocks in your protein:**")
    st.bar_chart(c["counts"])
    with st.expander("What does this mean?"):
        frac = c.get("class_fractions", {})
        top = max(frac, key=frac.get) if frac else "-"
        direction = ("stick to negatively charged surfaces"
                     if (ch or 0) > 0 else
                     "stick to positively charged surfaces"
                     if (ch or 0) < 0 else "stay neutral in solution")
        # Determine hydrophobicity label from GRAVY
    gravy_label = "hydrophilic" if gravy < 0 else "hydrophobic"

    msg = (
        "Your protein is mostly %s building blocks "
        "(%d in total, weighing about %d Da). "
        "It will tend to %s. "
        "The GRAVY (%.2f) indicates %s character, "
        "while the chart shows the amino-acid class distribution."
        ) % (top, c['length'], c['molecular_weight_da'], direction, gravy, gravy_label)

    st.write(msg)

    # NEW: Amino acid breakdown table with download
    from collections import Counter
    counts = Counter(seq)
    total = len(seq)
    sorted_aa = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    aa_df = {
        "Amino Acid": [a[0] for a in sorted_aa],
        "Count": [a[1] for a in sorted_aa],
        "Percentage": [f"{(a[1]/total*100):.1f}%" for a in sorted_aa]
    }
    st.dataframe(aa_df, hide_index=True, use_container_width=True)
    csv_lines = "Amino Acid,Count,Percentage\n"
    for aa, cnt in sorted_aa:
        csv_lines += f"{aa},{cnt},{cnt/total*100:.1f}\n"
    st.download_button(
        "⬇ Download AA table as CSV",
        csv_lines,
        file_name="amino_acid_counts.csv",
        mime="text/csv",
        help="Download a CSV with amino acid counts and percentages"
    )

