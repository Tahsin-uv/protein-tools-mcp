"""Structure tab for Protein Lab: upload PDB/mmCIF -> verdict cards."""
from __future__ import annotations

import os

import cif_tools
import pdb_utils

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")

_METAL_WORDS = {
    "ZN": "zinc", "FE": "iron", "MN": "manganese", "MG": "magnesium",
    "CU": "copper", "NI": "nickel", "CO": "cobalt", "CD": "cadmium",
}

_GEO_WORDS = {
    "linear": "a straight line of 2 atoms",
    "trigonal_planar": "a flat triangle of 3 atoms",
    "tetrahedral": "a pyramid of 4 atoms",
    "square_planar": "a flat square of 4 atoms",
    "trigonal_bipyramidal": "a double-pyramid of 5 atoms",
    "square_pyramidal": "a pyramid with a square base (5 atoms)",
    "octahedral": "an 8-faced shape of 6 atoms",
}


def _load(text: str, name: str):
    if name.lower().endswith(".cif") or text.lstrip().startswith("data_"):
        return cif_tools.to_pdb_atoms(cif_tools.parse_cif(text))
    return text


def _demo(name: str) -> str:
    with open(os.path.join(FIX, name), encoding="utf-8") as f:
        return f.read()

def render() -> None:
    import streamlit as st

    st.header("Is my structure healthy?")
    st.write("Upload a `.pdb` or `.cif` file — or explore the demo "
             "structures below.")

    up = st.file_uploader("Upload a PDB or mmCIF file",
                          type=["pdb", "cif", "ent"])
    if "_lab_struct" not in st.session_state:
        st.session_state["_lab_struct"] = None
    col1, col2 = st.columns(2)
    if col1.button("Demo: zinc protein (12 aa)"):
        st.session_state["_lab_struct"] = (_demo("tiny_zn.pdb"),
                                            "tiny_zn.pdb")
        st.rerun()
    if col2.button("Demo: iron protein (130 aa)"):
        st.session_state["_lab_struct"] = (
            _demo("A8ABE7_holo_Fe2.pdb"), "A8ABE7_holo_Fe2.pdb")
        st.rerun()

    text, name = None, None
    if up is not None:
        text = up.read().decode("utf-8", errors="replace")
        name = up.name
        st.session_state["_lab_struct"] = (text, name)
    elif st.session_state.get("_lab_struct"):
        text, name = st.session_state["_lab_struct"]
    if text is None:
        st.info("Upload a file or pick a demo structure to begin.")
        return

    try:
        struct = _load(text, name or "")
    except ValueError as exc:
        st.error(f"Couldn't parse that file: {exc}")
        return
    v = pdb_utils.validate_pdb(struct)
    st.subheader(f"Verdict for `{name}`")
    problems = list(v.get("format_errors", []) or [])
    problems += list(v.get("unknown_residues", []) or [])
    problems += list(v.get("non_integer_resseq", []) or [])
    counts = v.get("counts", {})
    if not problems:
        st.success(f"Your structure looks good — "
                   f"{counts.get('residues', '?')} building blocks in "
                   f"{len(counts.get('chains', []))} chain(s), no errors.")
    else:
        st.warning(f"Your structure parsed with "
                   f"{len(problems)} issue(s) — see details below.")

    m1, m2, m3 = st.columns(3)
    m1.metric("Building blocks", counts.get("residues", "?"))
    m2.metric("Atoms", counts.get("atoms", "?"))
    chains = counts.get("chains", [])
    m3.metric("Chains", ", ".join(chains) if chains else "?")

    gaps = v.get("chain_gaps", [])
    if gaps:
        st.warning(f"There are {len(gaps)} gap(s) in the chain — "
                   f"missing pieces the experiment couldn't see, e.g. "
                   + ", ".join(f"chain {g['chain']}: "
                                f"{g['missing']} missing after "
                                f"position {g['from']}" for g in gaps[:3])
                   + ("" if len(gaps) <= 3 else " ..."))
    else:
        st.caption("No chain gaps — the chain is continuous.")

    seqs = pdb_utils.sequences(struct)
    chains_seq = {k: d for k, d in seqs.items()
                  if not k.startswith("_")}
    if chains_seq:
        show = next(iter(chains_seq.values()))["sequence"]
        st.write("**First chain, one-letter code:**")
        st.code(show[:120] + ("..." if len(show) > 120 else ""))
        st.caption(f"Full length: {len(show)} letters. "
                   "Each letter is one amino acid.")

    st.subheader("Metal ions (helper atoms)")
    metals = _find_metals(struct)
    if not metals:
        st.write("No metal ions found — most student structures have none, "
                 "which is completely normal.")
    else:
        for el in metals:
            word = _METAL_WORDS.get(el, el.lower())
            with st.expander(f"{word.capitalize()} ({el}): "
                             f"{len(metals[el])} site(s)"):
                g = pdb_utils.metal_geometry(struct, el)
                for site in g.get("metal_sites", []):
                    _render_site(site, word)
    with st.expander("What does this mean?"):
        st.write("Some proteins hold small charged metal atoms (iron, "
                 "zinc, ...) that help them work — like a wrench held "
                 "in a hand. The 'grip' is the set of nearby protein "
                 "atoms touching the metal. A neat, symmetric grip "
                 "(e.g. a pyramid of 4) usually means a real binding "
                 "site; a loose or lopsided one deserves a second look "
                 "before you trust it.")


def _find_metals(struct) -> dict[str, list]:
    try:
        atoms = struct["atoms"] if isinstance(struct, dict) \
            else pdb_utils.parse_pdb(struct)["atoms"]
    except Exception:
        return {}
    out: dict[str, list] = {}
    for a in atoms:
        el = (a.get("element") or "").strip().upper()
        if el in pdb_utils.METAL_ELEMENTS:
            out.setdefault(el, []).append(a)
    return out


def _render_site(site: dict, word: str) -> None:
    import streamlit as st

    cn = site.get("coordination_number", 0)
    fit = (site.get("geometry_fit") or {})
    best = fit.get("best")
    shape = _GEO_WORDS.get(best, best or "an irregular arrangement")
    rmsd = fit.get("rmsd_deg")
    neat = rmsd is not None and rmsd < 15.0
    icon = "Healthy grip" if neat else "Check this grip"
    st.markdown(f"**{icon} — {site.get('site')}:** the {word} atom is "
                f"held by **{cn}** protein atoms in {shape}"
                + (f" (fits the ideal shape within {rmsd} degrees)."
                    if rmsd is not None else "."))
    for c in site.get("ligand_contacts", [])[:6]:
        st.caption(f"- {c.get('residue')}, atom {c.get('atom')} "
                   f"({c.get('dist_A')} A away, {c.get('type')})")
    if site.get("bidentate_chelates"):
        st.caption("Two of these contacts come from the same acidic "
                   "building block gripping with both oxygens — a "
                   "strong, common hold.")

