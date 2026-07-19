# -*- coding: utf-8 -*-
"""
Master script for the corrected PHB C1-C3 Zenodo reproducibility package.

Default mode is fast and validates the packaged CSV data against the numerical
claims used in the article. Full mode recomputes the heavy reduced scalar
Helmholtz diagnostics.

Examples
--------
Fast validation:
    python scripts/PHB_C1C3_FINAL_ZENODO_MASTER_SCRIPT.py

Full recomputation (can take tens of minutes):
    python scripts/PHB_C1C3_FINAL_ZENODO_MASTER_SCRIPT.py --full

Smoke tests for additional scripts:
    python scripts/PHB_C1C3_FINAL_ZENODO_MASTER_SCRIPT.py --smoke

Model status: reduced axisymmetric scalar Helmholtz only. Not full Maxwell.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
REF = ROOT / 'reference_results'
REFFIG = ROOT / 'reference_figures'

REQUIRED_REFERENCE = [
    REF / 'primary_mode_grid_extension.csv',
    REF / 'C3_deduplicated_distinct_modes_70x48.csv',
    REF / 'horn_leakage_empirical_grid_variability_floor.csv',
    REF / 'geometry_matching_metrics.csv',
    REF / 'table8_representative_type_size_family.csv',
    REFFIG / 'fig_C2_eta100_grid_stability.png',
    REFFIG / 'fig_leakage_vs_floor_extended.png',
]

def run(cmd):
    print('\nRUN:', ' '.join(map(str, cmd)), flush=True)
    subprocess.run([sys.executable, *map(str, cmd)], cwd=str(ROOT), check=True)

def fast_validate():
    missing = [str(p.relative_to(ROOT)) for p in REQUIRED_REFERENCE if not p.exists()]
    if missing:
        print('Missing required reference files:')
        for m in missing: print('  ', m)
        raise SystemExit(1)
    run([SCRIPTS / 'verify_article_numbers.py'])

def full_recompute():
    run([SCRIPTS / 'PHB_C1C3_final_strict_extension_calc.py'])
    run([SCRIPTS / 'finish_final_strict_extension.py'])
    print('\nFull recomputation finished. Outputs are in run_full_strict_extension/ and recomputed_results/.')

def smoke_tests():
    # Lightweight smoke mode: syntax-check every script and run a small numerical core search.
    # The stress-test and continuous-scan scripts are intentionally not run here, because even
    # their reduced numerical modes can be slow on weak CPUs. Run them explicitly when needed.
    import py_compile
    for script in sorted(SCRIPTS.glob('*.py')):
        print('PY_COMPILE:', script.name)
        py_compile.compile(str(script), doraise=True)
    run([SCRIPTS / 'PHB_C2C3_hyperbolic_uniqueness_search.py', '--preset', 'quick', '--outdir', ROOT / 'smoke_core', '--package', ROOT / 'smoke_core.zip'])
    print('\nSmoke tests completed: all scripts compiled and the core solver quick run passed.')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--full', action='store_true', help='recompute the heavy scalar Helmholtz calculations')
    ap.add_argument('--smoke', action='store_true', help='run lightweight smoke tests for all standalone scripts')
    args = ap.parse_args()
    fast_validate()
    if args.smoke:
        smoke_tests()
    if args.full:
        full_recompute()
    print('\nMASTER RESULT: PASS')

if __name__ == '__main__':
    main()
