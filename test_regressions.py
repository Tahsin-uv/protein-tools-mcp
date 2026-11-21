'''Regression tests for review findings 1-4 (reviewed commit ebfd573).

Offline (no API key, no network).  Run: python test_regressions.py [name-filter]
  R1 entry points: uvx / setup_client must name an executable that [project.scripts] provides.
  R2 pdb_by_id:    real PDB IDs (they contain digits) are accepted; local files are found
                   outside the install directory.
  R3 metal_geometry: reference angle lists are ascending; perfect CN 4-8 sites fit with ~0 RMSD.
  R4 doc drift:    keyword args shown in SKILL.md exist; asking for the wrong metal is visible.
'''
from __future__ import annotations

import asyncio
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
for _var in ('PROTEIN_TOOLS_API_KEY', 'PROTEIN_TOOLS_RELAY_URL', 'PROTEIN_TOOLS_STRUCTURE_DIR'):
    os.environ.pop(_var, None)
os.environ['PROTEIN_TOOLS_HOME'] = tempfile.mkdtemp(prefix='ptm_regress_')

import pdb_utils  # noqa: E402
import server  # noqa: E402

NL, DQ = chr(10), chr(34)

def _het(serial, name, resn, resseq, xyz, el):
    x, y, z = xyz
    return (f'HETATM{serial:5d} {name:<4s} {resn:>3s} A{resseq:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {el:>2s}')

def _site_pdb(metal, vectors, radius):
    lines = [_het(1, metal, metal, 1, (0.0, 0.0, 0.0), metal)]
    for k, v in enumerate(vectors, start=2):
        n = math.sqrt(sum(c * c for c in v))
        lines.append(_het(k, 'O', 'HOH', k, tuple(radius * c / n for c in v), 'O'))
    return NL.join(lines) + NL + 'END' + NL

def _write(text, name='site.pdb', folder=None):
    folder = Path(folder) if folder else Path(tempfile.mkdtemp(prefix='ptm_site_'))
    path = folder / name
    path.write_text(text, encoding='utf-8')
    return path

def _ring(n):
    return [(math.cos(2 * math.pi * k / n), math.sin(2 * math.pi * k / n), 0.0) for k in range(n)]

_AXIAL = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)]
_OCT = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
_TET = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
_CUBE = [(a, b, c) for a in (1, -1) for b in (1, -1) for c in (1, -1)]
IDEAL_SITES = {
    'tetrahedral': ('ZN', _TET, 2.39, 4),
    'square_planar': ('ZN', _OCT[:4], 2.39, 4),
    'trigonal_bipyramidal': ('FE', _AXIAL + _ring(3), 2.39, 5),
    'square_pyramidal': ('FE', _OCT[:4] + [(0, 0, 1)], 2.39, 5),
    'octahedral': ('FE', _OCT, 2.39, 6),
    'pentagonal_bipyramidal': ('CA', _AXIAL + _ring(5), 2.39, 7),
    'cubic': ('CA', _CUBE, 2.39, 8),
    'hexagonal_bipyramidal': ('CA', _AXIAL + _ring(6), 2.39, 8),
}
# --------------------------------------------------------------------- R1
def _scripts():
    out, in_table = {}, False
    for raw in (HERE / 'pyproject.toml').read_text(encoding='utf-8').splitlines():
        line = raw.split('#', 1)[0].strip()
        if line.startswith('['):
            in_table = line == '[project.scripts]'
        elif in_table and '=' in line:
            key, val = line.split('=', 1)
            out[key.strip().strip(DQ)] = val.strip().strip(DQ)
    assert out, 'no [project.scripts] table in pyproject.toml'
    return out

def test_r1_setup_client_registers_a_real_executable():
    scripts = _scripts()
    assert scripts.get('protein-tools-mcp') == 'server:main', \
        f'need a protein-tools-mcp script pointing at server:main; have {sorted(scripts)}'
    assert scripts.get('protein-tools') == 'server:main', 'keep the legacy protein-tools script'
    import setup_client
    entry = setup_client._entry(None)
    args = list(entry['args'])
    if entry['command'] == 'uvx':
        exe = args[args.index('--from') + 2] if '--from' in args else args[0]
    else:
        exe = Path(entry['command']).name
    assert exe in scripts, f'setup_client registers {exe!r} but the package only provides {sorted(scripts)}'

def test_r1_docs_do_not_name_a_missing_executable():
    scripts, bad = _scripts(), []
    for doc in ('README.md', 'pyproject.toml'):
        for line in (HERE / doc).read_text(encoding='utf-8').splitlines():
            words = line.lstrip('# ').split()
            if len(words) > 1 and words[0] == 'uvx' and not words[1].startswith('-') \
                    and words[1] not in scripts:
                bad.append(f'{doc}: {line.strip()}')
    assert not bad, '; '.join(bad)

# --------------------------------------------------------------------- R2
def test_r2_pdb_id_validation_accepts_real_ids():
    for pid in ['1ABC', '4HHB', '1abc', ' 4hhb ', 'ABCD', 'A8AB', '12AB']:
        err = json.loads(server.pdb_by_id(pid)).get('error')
        assert err is None, f'{pid!r} wrongly rejected: {err}'
    non_ascii = [chr(0xC9) + '1AB', ''.join(chr(0x661 + i) for i in range(4))]
    for pid in ['', 'ABC', 'ABCDE', '1AB!', 'AB C', '../ab'] + non_ascii:
        err = json.loads(server.pdb_by_id(pid)).get('error', '')
        assert 'pdb_id must be' in err, f'{pid!r} wrongly accepted'

def test_r2_local_lookup_outside_install_dir():
    folder = Path(tempfile.mkdtemp(prefix='ptm_pdbid_'))
    _write(_site_pdb('ZN', _TET, 2.0), 'qxyz.pdb', folder)
    os.environ['PROTEIN_TOOLS_STRUCTURE_DIR'] = str(folder)
    try:
        res = json.loads(server.pdb_by_id('QXYZ'))
    finally:
        os.environ.pop('PROTEIN_TOOLS_STRUCTURE_DIR', None)
    assert res['local']['present'] is True, f'PROTEIN_TOOLS_STRUCTURE_DIR ignored: {res}'
    old = os.getcwd()
    os.chdir(folder)
    try:
        res = json.loads(server.pdb_by_id('qxyz'))
    finally:
        os.chdir(old)
    assert res['local']['present'] is True, f'current directory not searched: {res}'
    res = json.loads(server.pdb_by_id('QXYZ'))
    assert res['local']['present'] is False, f'unexpected hit: {res}'

# --------------------------------------------------------------------- R3
def test_r3_reference_angle_lists_are_ascending():
    for cn in range(2, 9):
        for name, ideal in pdb_utils._ideal_geometries(cn).items():
            assert len(ideal) == cn * (cn - 1) // 2, f'CN{cn} {name}: {len(ideal)} angles'
            assert ideal == sorted(ideal), \
                f'CN{cn} {name}: reference angles must be ascending (observed angles are sorted)'

def test_r3_perfect_sites_fit_their_own_geometry():
    bad = []
    for name, (metal, vectors, radius, cn) in IDEAL_SITES.items():
        path = _write(_site_pdb(metal, vectors, radius), name + '.pdb')
        site = json.loads(server.metal_geometry(str(path), metal_element=metal, cutoff_angstrom=2.4))['metal_sites'][0]
        fit = site['geometry_fit']
        if site['coordination_number'] != cn or fit['best'] != name or fit['rmsd_deg'] > 0.5:
            bad.append(f'{name}: CN={site["coordination_number"]} best={fit["best"]} rmsd={fit["rmsd_deg"]}')
    assert not bad, '; '.join(bad)

# --------------------------------------------------------------------- R4
def test_r4_docs_use_real_keyword_arguments():
    tools = {t.name: set(t.inputSchema.get('properties', {}))
             for t in asyncio.run(server.mcp.list_tools())}
    bad = []
    for doc in ('SKILL.md', 'README.md', 'LISTING.md'):
        path = HERE / doc
        if not path.exists():
            continue
        for m in re.finditer('`([a-z_]+)[(]([^`]*)[)]`', path.read_text(encoding='utf-8')):
            name, argtxt = m.group(1), m.group(2)
            if name not in tools:
                continue
            for kw in re.findall('(?:^|[, (])([a-z_]+)[ ]*=', argtxt):
                if kw not in tools[name]:
                    bad.append(f'{doc}: {name}({kw}=...) but parameters are {sorted(tools[name])}')
    assert not bad, '; '.join(bad)

def test_r4_wrong_metal_is_visible_in_the_output():
    path = _write(_site_pdb('FE', _OCT, 2.1))
    res = json.loads(server.metal_geometry(str(path), metal_element='ZN'))
    assert res['n_sites'] == 0, res
    assert res.get('metals_present') == {'FE': 1}, f'metals_present missing or wrong: {res}'
    assert 'FE' in res.get('note', ''), f'note should name the metals that ARE present: {res}'
    res = json.loads(server.metal_geometry(str(path), metal_element='FE'))
    assert res['n_sites'] == 1 and res['metals_present'] == {'FE': 1}, res

TESTS = [
    test_r1_setup_client_registers_a_real_executable,
    test_r1_docs_do_not_name_a_missing_executable,
    test_r2_pdb_id_validation_accepts_real_ids,
    test_r2_local_lookup_outside_install_dir,
    test_r3_reference_angle_lists_are_ascending,
    test_r3_perfect_sites_fit_their_own_geometry,
    test_r4_docs_use_real_keyword_arguments,
    test_r4_wrong_metal_is_visible_in_the_output,
]

def main() -> int:
    only = sys.argv[1] if len(sys.argv) > 1 else ''
    failed = 0
    for fn in TESTS:
        if only and only not in fn.__name__:
            continue
        try:
            fn()
            print('PASS ', fn.__name__)
        except Exception as exc:
            failed += 1
            print('FAIL ', fn.__name__)
            print('     ', type(exc).__name__ + ': ' + str(exc)[:500])
    print()
    print('ALL REGRESSION TESTS PASSED' if not failed else f'{failed} REGRESSION TEST(S) FAILED')
    return 1 if failed else 0

if __name__ == '__main__':
    sys.exit(main())