"""Hugging Face Spaces entry point — thin shim over protein_lab.py.

Spaces runs `streamlit run app.py`; local dev runs
`streamlit run protein_lab.py` (or the `protein-lab` console script).
Both share the same code.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import protein_lab  # noqa: F401  (runs the app on import)
