"""Data-only sub-package holding the mini GROMACS run fixture.

Exists only so setuptools treats `fixtures/mdrun_apo` as a package and ships
its files (md.log, rmsd.xvg, state.cpt) next to the other demo data.
"""
