"""Data-only package holding the Protein Lab demo fixtures.

This is NOT an importable API - the directory exists so that setuptools
installs it next to the flat top-level modules. `lab_seq`, `lab_structure`
and `gen_fixtures` resolve demo data as ``<module dir>/fixtures``, which
must keep working after a pip install (site-packages/fixtures/).
"""
