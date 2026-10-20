"""`protein-lab` console entry point: launch the Protein Lab web UI locally."""
from __future__ import annotations

import os
import sys


def main() -> None:
    """Run `streamlit run protein_lab.py` from the installed package."""
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print("Protein Lab needs streamlit: pip install streamlit")
        raise SystemExit(1)
    app = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "protein_lab.py")
    sys.argv = ["streamlit", "run", app, *sys.argv[1:]]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
