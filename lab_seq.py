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

    st.header("What is my sequence?")
    st.write("Paste letters (like `MKLIVFWYHEPQDNRS`), a FASTA record "
             "(starting with `>`), or upload a `.fasta` file.")

    if "_lab_seq" not in st.session_state:
        st.session_state["_lab_seq"] = ""
    src = st.text_area("Sequence or FASTA", height=120,
                       placeholder=">my_protein\nMKLIVFWYHEPQDNRSACDE",
                       key="_lab_seq")
    up = st.file_uploader("...or upload a .fasta file",
                          type=["fasta", "fa", "txt"])
    if up is not None:
        src = up.read().decode("utf-8", errors="replace")
        st.session_state["_lab_seq"] = src
    if st.button("Try demo data instead"):
        st.session_state["_lab_seq"] = _demo_text()
        st.rerun()
    src = st.session_state.get("_lab_seq", "")
    if not (src or "").strip():
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

    with st.expander("Glossary - what do these words mean?"):
        for term, meaning in _GLOSSARY.items():
            st.markdown(f"**{term}:** {meaning}")


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
    c1.metric("Length", f"{c['length']} amino acids")
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
    st.write("**Building blocks in your protein:**")
    st.bar_chart(c["counts"])
    with st.expander("What does this mean?"):
        frac = c.get("class_fractions", {})
        top = max(frac, key=frac.get) if frac else "-"
        direction = ("stick to negatively charged surfaces"
                     if (ch or 0) > 0 else
                     "stick to positively charged surfaces"
                     if (ch or 0) < 0 else "stay neutral in solution")
        st.write(
            f"Your protein is mostly **{top}** building blocks "
            f"({c['length']} in total, weighing about "
            f"{c['molecular_weight_da']:,.0f} Da). "
            f"It will tend to {direction}. "
            f"The chart shows how many of each of the 20 amino-acid "
            f"types appear.")

