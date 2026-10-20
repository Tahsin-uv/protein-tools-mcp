# Third-party data notices

protein-tools consumes public scientific databases. Attribution and terms:

| Source | Data used | License / terms |
|---|---|---|
| RCSB PDB | coordinate files (mmCIF), entry metadata | PDB data is released into the public domain (CC0 1.0). Courtesy citation of the wwPDB is expected. |
| UniProtKB | protein entries, sequences, annotations | CC BY 4.0 — **attribution required**. Responses are cached (24 h) by the optional hosted relay, so redistributed or displayed data must credit UniProt. |
| AlphaFold DB | predicted structures and pLDDT | CC BY 4.0 — attribution required. |
| PROPKA | parsed `.pka` output only (no code bundled) | see the PROPKA project license. |

Notes:
- The PDB-wide metal-site index is a *derived* dataset: coordinates are CC0,
  but the derived coordination/geometry annotations are original work of the
  project author and are not public-domain.
- Cached upstream responses must never be presented as the author's own data.

Suggested citations:
- wwPDB / RCSB PDB: Berman et al., *Acta Cryst.* D66, 12–21 (2010).
- UniProt Consortium, *Nucleic Acids Res.* (current release).
- Varadi et al., *Nucleic Acids Res.* 50, D439–D444 (2022) (AlphaFold DB).
