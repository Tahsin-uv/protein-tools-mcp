"""Sequence-type detection and FASTA parsing.

Detects protein / DNA / RNA / ambiguous / unknown from alphabet composition,
with explicit T-vs-U logic and confidence notes rather than silent guesses.

Stdlib-only (Biopython optional; not required).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

# ---------------------------------------------------------------------------
# IUPAC alphabets
# ---------------------------------------------------------------------------

AA_ALPHABET: set[str] = frozenset(
    "ACDEFGHIKLMNPQRSTVWY" "BZJ" "XU" "O" "*-"
)

DNA_ALPHABET: set[str] = frozenset("ACGT" "RYSWKMBDHVN")
RNA_ALPHABET: set[str] = frozenset("ACGU" "RYSWKMBDHVN")
NUCLEIC_ALPHABET: set[str] = DNA_ALPHABET | RNA_ALPHABET
ALL_ALPHABET: set[str] = AA_ALPHABET | NUCLEIC_ALPHABET

TypeLabel = Literal["protein", "dna", "rna", "ambiguous", "unknown"]


@dataclass
class SequenceTypeDetail:
    has_T: bool = False
    has_U: bool = False
    protein_score: int = 0
    dna_score: int = 0
    rna_score: int = 0
    nonstandard: list[str] = field(default_factory=list)
    ambiguous_pair: bool = False
    note: str = ""
    confidence: Literal["high", "medium", "low"] = "low"


def classify_sequence_type(seq: str) -> tuple[TypeLabel, SequenceTypeDetail]:
    """Detect protein / DNA / RNA / ambiguous / unknown from alphabet.

    Returns (label, detail) so callers get both the label and the reasoning.
    """
    s = seq.strip()
    if not s:
        d = SequenceTypeDetail()
        d.note = "Empty sequence."
        d.confidence = "low"
        return ("unknown", d)

    upper = s.upper()
    has_T = "T" in upper
    has_U = "U" in upper
    detail: SequenceTypeDetail = SequenceTypeDetail(has_T=has_T, has_U=has_U)

    # Only treat T/U as nucleic-acid markers when the sequence is actually
    # compatible with a nucleic-acid alphabet. A protein-like string that
    # happens to contain T (Thr) or U (selenocysteine) must not be claimed
    # as DNA/RNA: e.g. "MKLIVFWYHEPQDNRS" contains T but also E/F/H/I/L/M/P/Q/W.
    nucleic_compatible = not (set(upper) - NUCLEIC_ALPHABET)
    if has_T and not has_U and nucleic_compatible:
        d = _score(upper, detail)
        d.confidence = "high"
        d.note = "T present, no U: classified as DNA."
        return ("dna", d)
    if has_U and not has_T and nucleic_compatible:
        d = _score(upper, detail)
        d.confidence = "high"
        d.note = "U present, no T: classified as RNA."
        return ("rna", d)
    if has_T and has_U:
        detail.ambiguous_pair = True
        _score(upper, detail)
        detail.confidence = "high"
        detail.note = "Both T and U present; cannot cleanly determine DNA vs RNA."
        return ("ambiguous", detail)

    d = _score(upper, detail)
    d.confidence = "medium"
    if d.dna_score == d.rna_score and d.dna_score == 0:
        if _is_protein_like(upper, detail):
            d.confidence = "high"
            d.note = "No T/U; non-nucleic-acid characters present; classified as protein."
            return ("protein", d)
    if d.dna_score == d.rna_score and d.dna_score > 0:
        d.confidence = "low"
        d.note = (
            "No T or U present; equal DNA/RNA compatibility; "
            "defaulting to DNA (could be either)."
        )
        return ("dna", d)
    if d.dna_score > d.rna_score:
        d.confidence = "medium"
        d.note = "No T/U; DNA-alphabet match exceeds RNA; classified as DNA."
        return ("dna", d)
    if d.rna_score > d.dna_score:
        d.confidence = "medium"
        d.note = "No T/U; RNA-alphabet match exceeds DNA; classified as RNA."
        return ("rna", d)

    if _is_protein_like(upper, detail):
        d.confidence = "high"
        d.note = (
            "No T/U; non-nucleic-acid characters present; "
            "classified as protein."
        )
        return ("protein", d)

    if detail.nonstandard:
        d.confidence = "low"
        d.note = "Contains characters outside protein and nucleic-acid alphabets."
        return ("unknown", d)

    d.note = "Could not classify."
    return ("unknown", d)


def _score(upper: str, detail: SequenceTypeDetail) -> SequenceTypeDetail:
    prot = dna = rna = 0
    nonstandard: list[str] = []
    for ch in upper:
        if ch in AA_ALPHABET:
            prot += 1
        elif ch in DNA_ALPHABET:
            dna += 1
        elif ch in RNA_ALPHABET:
            rna += 1
        else:
            nonstandard.append(ch)
    detail.protein_score = prot
    detail.dna_score = dna
    detail.rna_score = rna
    detail.nonstandard = list(dict.fromkeys(nonstandard))
    return detail


def _is_protein_like(upper: str, detail: SequenceTypeDetail) -> bool:
    protein_only = AA_ALPHABET - NUCLEIC_ALPHABET
    seen_any = any(ch in protein_only for ch in upper)
    return seen_any or (
        detail.protein_score > detail.dna_score
        and detail.dna_score == detail.rna_score
    )


def classify_sequence(seq: str) -> dict:
    """Public-facing dict version of classify_sequence_type."""
    label, detail = classify_sequence_type(seq)
    return {
        "type": label,
        "kind": "protein" if label == "protein" else (
            "nucleotide" if label in ("dna", "rna") else "unknown"
        ),
        "length": len(seq.strip()),
        "first_few_chars": seq.strip()[:30],
        "detail": {
            "has_T": detail.has_T,
            "has_U": detail.has_U,
            "protein_score": detail.protein_score,
            "dna_score": detail.dna_score,
            "rna_score": detail.rna_score,
            "nonstandard_chars": detail.nonstandard,
            "ambiguous_pair": detail.ambiguous_pair,
            "note": detail.note,
            "confidence": detail.confidence,
        },
    }


def parse_fasta(source: str, is_file: bool = True) -> list[dict]:
    """Parse one or more FASTA records from a file path or raw text.

    Returns a list of records:
        {
          "id": "<header> or first whitespace-delimited token",
          "description": "<rest of header after the id>",
          "sequence": "<full uppercase sequence, no whitespace>",
          "length": <int>,
          "type": {...from classify_sequence...},
          "maybe_uniprot_accession": "<accession or None, heuristic>"
        }

    If the source is a file path and the file does not exist, raises
    FileNotFoundError; if the content has no FASTA records, raises ValueError.
    """
    if is_file:
        try:
            with open(source, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"FASTA file not found: {source}")
        except OSError as e:
            raise OSError(f"Error reading FASTA file {source}: {e}")
    else:
        text = source

    records = _parse_fasta_records(text)
    if not records:
        raise ValueError("No FASTA records found in the provided source.")

    out: list[dict] = []
    for r in records:
        seq = r["sequence"].upper()
        out.append({
            "id": r["id"],
            "description": r["description"],
            "sequence": seq,
            "length": len(seq),
            "type": classify_sequence(seq),
            "maybe_uniprot_accession": _hint_uniprot_accession(r["id"]),
        })
    return out


def _hint_uniprot_accession(header_id: str) -> str | None:
    """Heuristic: return a UniProt-like accession from a FASTA header id.

    Handles 'sp|P00533|NAME' / 'tr|...' forms and bare accessions
    (5-10 alphanumerics starting with a letter). Returns None otherwise.
    """
    h = (header_id or "").strip()
    if "|" in h:
        parts = h.split("|")
        # 'sp|P00533|NAME' -> P00533; 'A0A001|NAME' -> A0A001 (bare|desc).
        cand = parts[1] if len(parts) >= 3 else parts[0]
        if cand and cand[0].isalpha():
            return cand
        return None
    bare = h.split()[0] if h else ""
    if 5 <= len(bare) <= 10 and bare[0].isalpha() and bare.isalnum():
        return bare
    return None


def _parse_fasta_records(text: str) -> list[dict]:
    records: list[dict] = []
    current_id: str | None = None
    current_desc: str | None = None
    current_seq_parts: list[str] = []
    in_seq = False

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if in_seq and current_id is not None:
                _flush(current_id, current_desc, current_seq_parts, records)
            header = line[1:].strip()
            parts = header.split(None, 1)
            current_id = parts[0] if parts else "unknown"
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq_parts = []
            in_seq = True
        else:
            if in_seq:
                current_seq_parts.append(line)

    if in_seq and current_id is not None:
        _flush(current_id, current_desc, current_seq_parts, records)
    return records


def _flush(cur_id: str, cur_desc: str | None, parts: list[str], records: list[dict]) -> None:
    seq = "".join(parts).upper()
    records.append({
        "id": cur_id,
        "description": cur_desc or "",
        "sequence": seq,
    })


def normalize_sequence(seq: str, target_type: str | None = None) -> str:
    """Clean a sequence: uppercase, strip gaps/whitespace.

    target_type: if 'rna' and the sequence contains 'T', convert T->U;
    if 'dna' and the sequence contains 'U', convert U->T.
    """
    s = seq.upper()
    allowed_gap = {" ", "\t", "\n", "\r", "-", "*", ".", "/", ":", ";"}
    kept: list[str] = []
    for ch in s:
        if ch in allowed_gap:
            continue
        kept.append(ch)
    result = "".join(kept)
    if target_type == "rna" and "T" in result:
        result = result.replace("T", "U")
    elif target_type == "dna" and "U" in result:
        result = result.replace("U", "T")
    return result


# ---------------------------------------------------------------------------
# Protein physicochemistry (Tool 1: aa_composition)
# ---------------------------------------------------------------------------

# Average isotopic masses (Da); water added once for termini.
AA_MASS: dict[str, float] = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
    "U": 150.0388, "O": 255.3172,
}
_WATER_MASS = 18.01528

# Kyte-Doolittle hydrophobicity scale.
GRAVY_SCALE: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
    "K": -3.9, "R": -4.5,
}

# Physicochemical classes (Cys in its own thiol bin).
AA_CLASS: dict[str, str] = {
    "A": "hydrophobic", "V": "hydrophobic", "I": "hydrophobic",
    "L": "hydrophobic", "M": "hydrophobic", "F": "hydrophobic",
    "W": "hydrophobic", "Y": "aromatic", "P": "polar",
    "G": "polar", "S": "polar", "T": "polar", "N": "polar",
    "Q": "polar", "C": "cysteine", "K": "basic", "R": "basic",
    "H": "basic", "D": "acidic", "E": "acidic",
}

# pKa set for the pI estimate (Bjellqvist-derived, termini included).
_PKA = {"N_term": 8.0, "C_term": 3.1, "D": 3.9, "E": 4.3, "H": 6.5,
        "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5}


def aa_composition(seq: str) -> dict:
    """Physicochemical summary of a protein sequence.

    Rejects non-protein input with ValueError naming the detected type.
    pI is a charge-curve bisection estimate (documented as approximate).
    """
    s = normalize_sequence(seq)
    if not s:
        raise ValueError("Empty sequence: nothing to analyze.")
    label, detail = classify_sequence_type(s)
    if label != "protein":
        raise ValueError(
            f"Not a protein sequence (detected: {label}; {detail.note})")
    counts: dict[str, int] = {}
    for ch in s:
        counts[ch] = counts.get(ch, 0) + 1
    n = len(s)
    standard = [ch for ch in s if ch in AA_MASS]
    mw = sum(AA_MASS[ch] for ch in standard) + _WATER_MASS
    classes: dict[str, int] = {}
    for ch in standard:
        c = AA_CLASS.get(ch, "other")
        classes[c] = classes.get(c, 0) + 1
    gravy_vals = [GRAVY_SCALE[ch] for ch in standard if ch in GRAVY_SCALE]
    # Net charge at pH 7: K/R +1, D/E -1, H +0.5 (documented estimate).
    charge7 = (counts.get("K", 0) + counts.get("R", 0)
               - counts.get("D", 0) - counts.get("E", 0)
               + 0.5 * counts.get("H", 0))
    return {
        "length": n,
        "molecular_weight_da": round(mw, 2),
        "counts": counts,
        "percent": {k: round(100 * v / n, 1) for k, v in counts.items()},
        "class_counts": classes,
        "class_fractions": {k: round(v / len(standard), 3)
                            for k, v in classes.items()} if standard else {},
        "charge_at_pH7_estimate": round(charge7, 1),
        "gravy": round(sum(gravy_vals) / len(gravy_vals), 3)
        if gravy_vals else None,
        "pi_estimate": _estimate_pi(counts),
        "ambiguous_letters": sorted(set(s) - set(AA_MASS)),
        "note": ("pI is a charge-curve bisection estimate (Bjellqvist-like "
                 "pKa set); MW excludes ambiguous letters; charge at pH 7 "
                 "uses K/R +1, D/E -1, H +0.5."),
    }


def _estimate_pi(counts: dict[str, int]) -> float:
    """Isoelectric point by bisection of the net-charge curve (pH 0-14)."""

    def charge(ph: float) -> float:
        pos = (1.0 / (1.0 + 10 ** (ph - _PKA["N_term"])))
        pos += sum(counts.get(aa, 0) / (1.0 + 10 ** (ph - _PKA[aa]))
                   for aa in ("K", "R", "H"))
        neg = (1.0 / (1.0 + 10 ** (_PKA["C_term"] - ph)))
        neg += sum(counts.get(aa, 0) / (1.0 + 10 ** (_PKA[aa] - ph))
                   for aa in ("D", "E", "C", "Y"))
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


# ---------------------------------------------------------------------------
# Multi-FASTA summary (Tool 2: fasta_stats)
# ---------------------------------------------------------------------------

def fasta_stats(source: str, is_file: bool = True) -> dict:
    """Assembly-style summary over 1+ FASTA records.

    GC% is computed over DNA (+ambiguous) records only; null with a note
    when no DNA records exist (never 0-by-default).
    """
    records = parse_fasta(source, is_file=is_file)
    lengths = sorted(r["length"] for r in records)
    total = sum(lengths)
    # N50: shortest record such that records >= it cover half the residues.
    half, n50, run = total / 2, lengths[-1], 0
    for ln in sorted(lengths, reverse=True):
        run += ln
        if run >= half:
            n50 = ln
            break
    by_type: dict[str, int] = {}
    for r in records:
        t = r["type"]["type"]
        by_type[t] = by_type.get(t, 0) + 1
    gc_rows = []
    for r in records:
        if r["type"]["type"] in ("dna", "ambiguous"):
            sq = r["sequence"]
            gc = sum(1 for ch in sq if ch in "GC")
            atgc = sum(1 for ch in sq if ch in "ACGTU")
            if atgc:
                gc_rows.append({"id": r["id"],
                               "gc_percent": round(100 * gc / atgc, 1)})
    n_gc = len(gc_rows)
    mean_gc = (round(sum(g["gc_percent"] for g in gc_rows) / n_gc, 1)
               if n_gc else None)
    nl = len(lengths)
    if nl % 2:
        median = float(lengths[nl // 2])
    else:
        median = round((lengths[nl // 2 - 1] + lengths[nl // 2]) / 2, 1)
    return {
        "n_records": len(records),
        "length_min": lengths[0],
        "length_max": lengths[-1],
        "length_mean": round(total / len(records), 1),
        "length_median": median,
        "n50": n50,
        "total_residues": total,
        "by_type": by_type,
        "gc_percent_mean": mean_gc,
        "n_gc_records": n_gc,
        "gc_note": (f"GC% over {n_gc} dna/ambiguous record(s)."
                    if n_gc else
                    "No DNA records: GC% not applicable (null, not 0)."),
        "gc_per_record": gc_rows,
        "records": [{"id": r["id"], "length": r["length"],
                     "type": r["type"]["type"]} for r in records],
    }
    """Heuristic: if the FASTA header looks like a UniProt accession,
    return it. Accepts sp|Pxxx|NAME, tr|Qxxx|NAME, A0Axxx, P12345, etc."""
    h = header_id.strip()
    for prefix in ("sp|", "tr|", "PE|", "UP|", "ID|"):
        if h.startswith(prefix):
            acc = h.split("|")[1] if "|" in h else ""
            if acc and acc[0].isalpha():
                return acc
    # bare accession: starts with letter, 5-10 alphanumerics
    bare = h.split("|")[0].split()[0]
    if 5 <= len(bare) <= 10 and bare[0].isalpha() and bare.isalnum():
        return bare
    return None
