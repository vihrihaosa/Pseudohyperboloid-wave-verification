# -*- coding: utf-8 -*-
"""
Fast numerical consistency checker for the PHB C1-C3 Zenodo package.

This script does not solve the eigenproblem. It checks the reference CSV tables
included in the Zenodo package against the numerical statements used in the
article. Use the master script with --full to recompute the heavy scalar
Helmholtz runs.

Model status: reduced axisymmetric scalar Helmholtz only; not full Maxwell.
"""
from __future__ import annotations
from pathlib import Path
import math
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'reference_results'

class CheckError(RuntimeError):
    pass

def _close(name: str, got: float, expected: float, tol: float) -> None:
    if not math.isfinite(got) or abs(got - expected) > tol:
        raise CheckError(f"{name}: got {got:.12g}, expected {expected:.12g} ± {tol}")
    print(f"OK  {name}: {got:.6g}")

def _pct(v):
    return 100.0 * float(v)

def check_primary_grid() -> None:
    df = pd.read_csv(RES / 'primary_mode_grid_extension.csv')
    targets = {
        ('PH', 70, 48): (15.7678, 93.026, 46.510, 23.892, 6.258, 5.519921e-05),
        ('PH', 100, 70): (15.8017, 92.501, 50.650, 25.054, 6.678, 2.098427e-05),
        ('PH', 130, 92): (15.8140, 92.444, 50.740, 25.346, 6.487, 3.254676e-05),
        ('PH', 160, 112): (15.8153, 92.380, 51.194, 25.820, 6.645, 2.264978e-05),
        ('LINEAR', 130, 92): (15.8870, 92.151, 50.186, 24.622, 4.948, 1.134427e-03),
        ('POLY2', 130, 92): (16.0522, 91.401, 49.812, 24.731, 2.263, 5.432239e-01),
        ('SMOOTHSTEP', 130, 92): (15.9811, 91.457, 49.909, 24.466, 3.719, 1.9105e-02),
    }
    for (geom, nx, ns), (kr, eta100, eta50, eta25, logcf, leak_pct) in targets.items():
        r = df[(df.geometry == geom) & (df.Nx == nx) & (df.Ns == ns)]
        if r.empty:
            raise CheckError(f"missing row {geom} {nx}x{ns}")
        r = r.iloc[0]
        _close(f"primary {geom} {nx}x{ns} kR", float(r.kR), kr, 5e-4)
        _close(f"primary {geom} {nx}x{ns} eta100%", float(r.eta_100R_pct), eta100, 0.015)
        _close(f"primary {geom} {nx}x{ns} eta50%", float(r.eta_050R_pct), eta50, 0.02)
        _close(f"primary {geom} {nx}x{ns} eta25%", float(r.eta_025R_pct), eta25, 0.02)
        _close(f"primary {geom} {nx}x{ns} log10CF", float(r.log10_CF), logcf, 0.004)
        _close(f"primary {geom} {nx}x{ns} leakage%", float(r.horn_leakage_pct), leak_pct, max(1e-8, leak_pct*0.02))

def check_c3() -> None:
    df = pd.read_csv(RES / 'C3_deduplicated_distinct_modes_70x48.csv')
    ph = df[df.geometry == 'PH'].sort_values('kR')
    above = ph[ph.eta_100R >= 0.70]
    if len(above) != 2:
        raise CheckError(f"PH C3 over-threshold distinct modes: got {len(above)}, expected 2")
    expected = [(15.7678, 93.026, 46.510), (17.6892, 81.454, 40.982)]
    for i, (kr, e100, e50) in enumerate(expected):
        r = above.iloc[i]
        _close(f"C3 PH peak {i+1} kR", float(r.kR), kr, 5e-4)
        _close(f"C3 PH peak {i+1} eta100%", _pct(r.eta_100R), e100, 0.02)
        _close(f"C3 PH peak {i+1} eta50%", _pct(r.eta_050R), e50, 0.02)
    # Also verify the two below-threshold rows that must not be silently omitted from a full table.
    for kr in [20.1143, 20.3628, 21.7768, 22.1449, 22.2592, 22.5921]:
        if (ph.kR.sub(kr).abs() < 8e-4).sum() == 0:
            raise CheckError(f"missing below-threshold PH C3 mode near kR={kr}")
    print('OK  C3: two over-threshold PH peaks and below-threshold intervening modes are present')

def check_leakage_floor() -> None:
    df = pd.read_csv(RES / 'horn_leakage_empirical_grid_variability_floor.csv')
    targets = {
        'CYLINDER_R': (23.9224, 1.008553, 23.7195),
        'ELLIPSOID_L': (1.328991, 0.309401, 4.29537),
        'HERMITE_SLOPE_MATCH': (5.125418e-08, 1.065397e-07, 0.481081),
        'LINEAR': (7.858354e-04, 5.235212e-04, 1.50106),
        'PH': (2.264978e-05, 3.421584e-05, 0.66197),
        'POLY2': (0.407577, 0.137978, 2.95392),
        'POLY_VOL_MATCH': (1.062060e-05, 1.662557e-05, 0.638817),
        'SMOOTHSTEP': (0.012660, 0.008699, 1.45530),
    }
    for geom, (fine, floor, ratio) in targets.items():
        r = df[df.geometry == geom]
        if r.empty:
            raise CheckError(f"missing leakage-floor row {geom}")
        r = r.iloc[0]
        _close(f"floor {geom} fine leakage%", float(r.fine_horn_leakage_pct), fine, max(1e-8, abs(fine)*0.002))
        _close(f"floor {geom} empirical floor%", float(r.empirical_grid_variability_floor_pct), floor, max(1e-8, abs(floor)*0.002))
        _close(f"floor {geom} ratio", float(r.fine_leakage_over_floor), ratio, max(1e-5, abs(ratio)*0.002))

def check_table8_family() -> None:
    # Table 8 values are a representative subset of the PH rows in primary_mode_grid_extension.csv.
    # In this corrected package they are explicitly stored as a small reference table.
    fp = RES / 'table8_representative_type_size_family.csv'
    df = pd.read_csv(fp)
    if len(df) != 5:
        raise CheckError('Table 8 reference family must contain 5 rows')
    if df['eta_100R_pct'].min() < 92.0:
        raise CheckError('Table 8 reference family contains eta_100R below 92%')
    print('OK  Table 8 reference family: 5 representative b-values above 92%')

def main() -> int:
    required = [
        'primary_mode_grid_extension.csv',
        'C3_deduplicated_distinct_modes_70x48.csv',
        'horn_leakage_empirical_grid_variability_floor.csv',
        'table8_representative_type_size_family.csv',
    ]
    for fn in required:
        if not (RES / fn).exists():
            raise CheckError(f"missing required reference result: {fn}")
    print('PHB C1-C3 Zenodo numerical-consistency check')
    print(f'Package root: {ROOT}')
    check_primary_grid()
    check_c3()
    check_leakage_floor()
    check_table8_family()
    print('\nRESULT: PASS. Reference CSV data match the corrected article-level numerical statements.')
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except CheckError as e:
        print('\nRESULT: FAIL')
        print(e)
        raise SystemExit(1)
