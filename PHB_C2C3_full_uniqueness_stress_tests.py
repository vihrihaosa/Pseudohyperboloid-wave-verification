# -*- coding: utf-8 -*-
"""
PHB_C2C3_full_uniqueness_stress_tests.py

Reduced scalar C2-C3 uniqueness stress tests for the second-order vertical
pseudohyperboloid (PHB).  This script extends the ordinary annular-energy
criterion eta(|x|<=c, |rho-R|<=0.10R) by tests aimed at the wave imprint of the
hyperbolic focal billiard:

  1. Hyperbolic-law perturbation test: does inter-focal confinement collapse
     when the horn is slightly deformed while keeping the same endpoints?
  2. Focal-plane scan: is the energy barrier tied to the true c=sqrt(a^2+b^2)?
  3. Cap-dependence test: does the effect survive changes of the central cap?
  4. m-sweep: does the signature grow with azimuthal order?
  5. scale-invariance smoke test: do dimensionless metrics survive scaling?

Model status:
  Boundary-fitted axisymmetric scalar Helmholtz eigenproblem.
  This is not full-vector Maxwell and not a laser calculation.

Recommended:
  python PHB_C2C3_full_uniqueness_stress_tests.py --preset standard

The preset 'quick' is a smoke test.  The preset 'standard' is the one used for
this report.  The preset 'confirm' adds a few heavier 100x70 confirmations.
"""
from __future__ import annotations

import argparse, importlib.util, json, math, time, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load the previous verified solver module so this script remains compact.
BASE_SCRIPT = Path(__file__).resolve().with_name('PHB_C2C3_hyperbolic_uniqueness_search.py')
spec = importlib.util.spec_from_file_location('phb_core', str(BASE_SCRIPT))
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

BASE_GEOMETRY_LABELS = core.geometry_labels
# -----------------------------------------------------------------------------
# Extended geometries: perturbed PHB horns and central-cap variants.
# -----------------------------------------------------------------------------

def _eps_from_kind(kind: str) -> float:
    # PH_PERT_PLUS_001 => 0.01, PH_PERT_MINUS_005 => 0.05, etc.
    code = kind.split('_')[-1]
    if code == '001': return 0.01
    if code == '002': return 0.02
    if code == '005': return 0.05
    if code == '010': return 0.10
    raise ValueError(kind)


def ext_horn(absx, a, b, R, L, kind):
    x = np.asarray(absx, dtype=float)
    t = np.clip((x - a) / max(L-a, 1e-14), 0.0, 1.0)
    # Underlying standard controls from the old script.
    if kind in ['PH','LINEAR','POLY0p5','POLY1p5','POLY2','POLY3','SMOOTHSTEP','CIRCULAR','COSINE']:
        return core._profile_horn_by_t(x, a, b, R, L, kind)
    # Central cap variants still use the exact hyperbolic horn.
    if kind.startswith('PH_CAP_'):
        return core._profile_horn_by_t(x, a, b, R, L, 'PH')
    # Perturb the PHB horn while preserving endpoints at x=a and x=L.
    if kind.startswith('PH_PERT_PLUS_') or kind.startswith('PH_PERT_MINUS_'):
        base = core._profile_horn_by_t(x, a, b, R, L, 'PH')
        eps = _eps_from_kind(kind)
        sign = 1.0 if 'PLUS' in kind else -1.0
        shape = np.sin(np.pi*t)**2  # zero at endpoints, strongest in horn middle
        return np.maximum(1e-8*R, base + sign*eps*R*shape)
    raise ValueError(f'Unknown geometry: {kind}')


def ext_wall_profile(xs, kind, a, b, R, L=None):
    if L is None:
        L = core.common_length(a,b,R)
    x = np.abs(np.asarray(xs, dtype=float))
    out = np.zeros_like(x)
    tor = x <= a
    # Central cap variants.  All are continuous at |x|=a and have max R+a at x=0
    # except FLAT_R, which deliberately removes the cap bulge.
    if kind == 'PH_CAP_FLAT_R':
        out[tor] = R
    elif kind == 'PH_CAP_PARABOLIC':
        out[tor] = R + a*(1.0 - (x[tor]/max(a,1e-14))**2)
    elif kind == 'PH_CAP_COSINE':
        out[tor] = R + a*0.5*(1.0 + np.cos(np.pi*x[tor]/max(a,1e-14)))
    else:
        out[tor] = R + np.sqrt(np.maximum(0.0, a*a - x[tor]**2))
    horn = (x > a) & (x <= L)
    out[horn] = ext_horn(x[horn], a,b,R,L,kind)
    return out


def ext_wall_derivative(xs, kind, a, b, R, L=None):
    if L is None:
        L = core.common_length(a,b,R)
    h = 1e-5*max(L,R,a,b)
    xs = np.asarray(xs, dtype=float)
    x1 = np.clip(xs-h, -L+1e-8, L-1e-8)
    x2 = np.clip(xs+h, -L+1e-8, L-1e-8)
    d = (ext_wall_profile(x2,kind,a,b,R,L)-ext_wall_profile(x1,kind,a,b,R,L))/np.maximum(x2-x1,1e-15)
    return np.clip(d, -300.0, 300.0)


def ext_geometry_labels():
    labels = BASE_GEOMETRY_LABELS()
    labels.update({
        'PH_PERT_PLUS_001':'PH horn +1% endpoint-preserving deformation',
        'PH_PERT_PLUS_002':'PH horn +2% endpoint-preserving deformation',
        'PH_PERT_PLUS_005':'PH horn +5% endpoint-preserving deformation',
        'PH_PERT_PLUS_010':'PH horn +10% endpoint-preserving deformation',
        'PH_PERT_MINUS_001':'PH horn -1% endpoint-preserving deformation',
        'PH_PERT_MINUS_002':'PH horn -2% endpoint-preserving deformation',
        'PH_PERT_MINUS_005':'PH horn -5% endpoint-preserving deformation',
        'PH_CAP_FLAT_R':'PH hyperbolic horn + flat R cap',
        'PH_CAP_PARABOLIC':'PH hyperbolic horn + parabolic cap',
        'PH_CAP_COSINE':'PH hyperbolic horn + cosine cap',
    })
    return labels

# Monkey-patch solver geometry functions.
core.wall_profile = ext_wall_profile
core.wall_derivative = ext_wall_derivative
core.geometry_labels = ext_geometry_labels

# -----------------------------------------------------------------------------
# Helper analysis utilities.
# -----------------------------------------------------------------------------

def best_mode_for_geometry(kind, a,b,R,m,Nx,Ns,kR_targets, bc='natural_TE_like_scalar_proxy', nev=6):
    env, df = core.solve_shift_scan(Nx,Ns,a,b,R,kind,bc,m,np.asarray(kR_targets,float),nev=nev,tol=1e-6,maxiter=2800)
    if df.empty:
        return None, env, df
    # For uniqueness we rank by the combined score, not by eta alone.
    best = df.sort_values('hyperbolic_signature_candidate_score', ascending=False).iloc[0].to_dict()
    return best, env, df


def get_axial_density(kind, a,b,R,m,Nx,Ns,kR_target, bc='natural_TE_like_scalar_proxy'):
    res = core.solve_single_target_with_density(Nx,Ns,a,b,R,kind,bc,m,kR_target,nev=8)
    if res is None:
        return None
    row, xs, dens, prob = res
    return row, xs, dens, prob


def c_factor_scan_from_density(row, xs, dens, a,b,R, factors):
    c = math.sqrt(a*a+b*b)
    dx = np.mean(np.diff(xs))
    out=[]
    for f in factors:
        cf = f*c
        inter = float(np.trapz(dens[np.abs(xs)<=cf], xs[np.abs(xs)<=cf])) if np.any(np.abs(xs)<=cf) else 0.0
        horn = max(1.0-inter, 1e-300)
        # strip-based drop using axial density only
        strip = max(0.05*cf, 2.5*dx)
        inside_mask = (np.abs(xs) >= cf-strip) & (np.abs(xs) <= cf)
        outside_mask = (np.abs(xs) > cf) & (np.abs(xs) <= cf+strip)
        Ein = float(np.trapz(dens[inside_mask], xs[inside_mask])) if np.any(inside_mask) else 0.0
        Eout = float(np.trapz(dens[outside_mask], xs[outside_mask])) if np.any(outside_mask) else 0.0
        out.append({
            'geometry':row['geometry'], 'm':row['m'], 'kR':row['kR'],
            'c_factor':f, 'c_used':cf, 'E_absx_le_cfactor_c':inter,
            'E_absx_gt_cfactor_c':horn, 'CF_cfactor':inter/max(horn,1e-300),
            'drop_ratio_cfactor':Ein/max(Eout,1e-300),
        })
    return pd.DataFrame(out)


def save_plot_bar(df, xcol, ycol, outpath, title, ylabel, logy=True, top=None):
    if df.empty: return
    d=df.copy()
    if top: d=d.head(top)
    plt.figure(figsize=(11,5.5))
    plt.bar(d[xcol].astype(str), d[ycol])
    if logy: plt.yscale('log')
    plt.xticks(rotation=40, ha='right')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def save_plot_line(df, xcol, ycols, outpath, title, ylabel, logy=False):
    if df.empty: return
    plt.figure(figsize=(9,5))
    for name, sub in df.groupby('geometry'):
        for ycol in ycols:
            plt.plot(sub[xcol], sub[ycol], marker='o', label=f'{name}: {ycol}')
    if logy: plt.yscale('log')
    plt.xlabel(xcol)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def run_all(args):
    t0=time.time()
    outdir=Path(args.outdir)
    resdir=outdir/'results'; figdir=outdir/'figures'
    resdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(parents=True, exist_ok=True)
    a,b,R=0.3,0.6,3.0
    bc='natural_TE_like_scalar_proxy'

    if args.preset == 'quick':
        # Very small smoke-test preset: verifies that the script and imports work.
        # It is not the standard stress-test run used for the article.
        Nx,Ns=28,20; m_main=15; targets=np.array([15.8])
        do_confirm=False
    elif args.preset == 'confirm':
        Nx,Ns=90,62; m_main=15; targets=np.arange(15.25,18.26,0.25)
        do_confirm=True
    else:
        Nx,Ns=70,48; m_main=15; targets=np.arange(15.25,18.26,0.25)
        do_confirm=True

    all_modes=[]; env_all=[]

    # Test 1: focal-law perturbation and non-hyperbolic controls.
    perturb_geoms=(['PH','LINEAR'] if args.preset=='quick' else ['PH','PH_PERT_PLUS_001','PH_PERT_PLUS_002','PH_PERT_PLUS_005','PH_PERT_PLUS_010',
                   'PH_PERT_MINUS_001','PH_PERT_MINUS_002','PH_PERT_MINUS_005','LINEAR','POLY2','SMOOTHSTEP'])
    perturb_rows=[]
    for kind in perturb_geoms:
        best, env, df = best_mode_for_geometry(kind,a,b,R,m_main,Nx,Ns,targets,bc=bc,nev=(3 if args.preset=='quick' else 5))
        if not env.empty: env_all.append(env)
        if not df.empty: all_modes.append(df)
        if best is not None:
            perturb_rows.append(best)
            print('PERT/CNTRL', kind, 'kR', f"{best['kR']:.4f}", 'eta', f"{100*best['eta_100R']:.2f}", 'CF', f"{best['CF_axial_inter_over_horns']:.3g}")
    df_pert=pd.DataFrame(perturb_rows).sort_values('CF_axial_inter_over_horns', ascending=False)
    df_pert.to_csv(resdir/'test1_perturbation_and_controls_best_modes.csv', index=False)

    # Test 2: focal-plane c-factor scan for main geometries using axial density.
    cfactors=np.arange(0.6,1.61,0.1)
    cscan=[]
    dens_rows=[]
    for kind in (['PH','LINEAR'] if args.preset=='quick' else ['PH','LINEAR','POLY2','SMOOTHSTEP']):
        # Use the best kR found above if available.
        if df_pert[df_pert.geometry==kind].empty:
            continue
        krt=float(df_pert[df_pert.geometry==kind].iloc[0]['kR'])
        res=get_axial_density(kind,a,b,R,m_main,Nx,Ns,krt,bc=bc)
        if res is None: continue
        row,xs,dens,prob=res
        dens_rows.append(pd.DataFrame({'geometry':kind,'x':xs,'axial_density':dens}))
        cscan.append(c_factor_scan_from_density(row,xs,dens,a,b,R,cfactors))
    df_cscan=pd.concat(cscan, ignore_index=True) if cscan else pd.DataFrame()
    df_cscan.to_csv(resdir/'test2_focal_plane_cfactor_scan.csv', index=False)
    if dens_rows:
        pd.concat(dens_rows, ignore_index=True).to_csv(resdir/'test2_axial_density_curves.csv', index=False)

    # Test 3: central cap dependence.
    cap_geoms=(['PH'] if args.preset=='quick' else ['PH','PH_CAP_FLAT_R','PH_CAP_PARABOLIC','PH_CAP_COSINE'])
    cap_rows=[]
    for kind in cap_geoms:
        best, env, df = best_mode_for_geometry(kind,a,b,R,m_main,Nx,Ns,targets,bc=bc,nev=(3 if args.preset=='quick' else 5))
        if not env.empty: env_all.append(env)
        if not df.empty: all_modes.append(df)
        if best is not None:
            cap_rows.append(best)
            print('CAP', kind, 'kR', f"{best['kR']:.4f}", 'eta', f"{100*best['eta_100R']:.2f}", 'CF', f"{best['CF_axial_inter_over_horns']:.3g}")
    df_cap=pd.DataFrame(cap_rows).sort_values('CF_axial_inter_over_horns', ascending=False)
    df_cap.to_csv(resdir/'test3_cap_dependence_best_modes.csv', index=False)

    # Test 4: m-sweep, PH vs LINEAR vs POLY2 vs SMOOTHSTEP.
    if args.preset=='quick':
        m_values=[15]
        Nx_m,Ns_m=28,20
    else:
        m_values=[6,10,12,15,20,25]
        Nx_m,Ns_m=60,42
    msweep=[]
    for m in m_values:
        mtargs=np.array([max(6.0,m*0.9), max(8.0,m*1.0), max(10.0,m*1.05), max(12.0,m*1.15)])
        for kind in (['PH','LINEAR'] if args.preset=='quick' else ['PH','LINEAR','POLY2','SMOOTHSTEP']):
            best, env, df=best_mode_for_geometry(kind,a,b,R,m,Nx_m,Ns_m,mtargs,bc=bc,nev=(3 if args.preset=='quick' else 4))
            if best is not None:
                msweep.append(best)
                print('M', m, kind, 'kR', f"{best['kR']:.3f}", 'eta', f"{100*best['eta_100R']:.1f}", 'CF', f"{best['CF_axial_inter_over_horns']:.2g}")
    df_ms=pd.DataFrame(msweep)
    if not df_ms.empty:
        df_ms.to_csv(resdir/'test4_m_sweep_best_modes.csv', index=False)

    # Test 5: scale invariance smoke test for PH.
    scales=([1.0] if args.preset=='quick' else [0.5,1.0,2.0])
    scale_rows=[]
    for s in scales:
        best, env, df=best_mode_for_geometry('PH',a*s,b*s,R*s,m_main,(28 if args.preset=='quick' else 55),(20 if args.preset=='quick' else 38),([15.8] if args.preset=='quick' else [15.5,15.8,16.0]),bc=bc,nev=(3 if args.preset=='quick' else 4))
        if best is not None:
            best['scale_factor']=s
            scale_rows.append(best)
    df_scale=pd.DataFrame(scale_rows)
    df_scale.to_csv(resdir/'test5_scale_invariance_smoke.csv', index=False)

    # Optional 100x70 confirmation for the strongest perturb/control/cap cases.
    confirm_rows=[]
    if do_confirm:
        confirm_geoms=['PH','LINEAR','POLY2','SMOOTHSTEP','PH_PERT_PLUS_005','PH_PERT_MINUS_005','PH_CAP_PARABOLIC','PH_CAP_FLAT_R']
        for kind in confirm_geoms:
            best, env, df=best_mode_for_geometry(kind,a,b,R,m_main,100,70,[15.55,15.80,16.05,17.80],bc=bc,nev=5)
            if best is not None:
                confirm_rows.append(best)
                print('CONFIRM100', kind, 'kR', f"{best['kR']:.4f}", 'eta', f"{100*best['eta_100R']:.2f}", 'CF', f"{best['CF_axial_inter_over_horns']:.3g}")
    df_confirm=pd.DataFrame(confirm_rows)
    if not df_confirm.empty:
        df_confirm.to_csv(resdir/'test6_100x70_selected_confirmation.csv', index=False)

    # Aggregate all modes and matched comparisons for standard geometries.
    if all_modes:
        df_all=pd.concat(all_modes, ignore_index=True)
        df_all.to_csv(resdir/'all_distinct_modes_from_stress_tests.csv', index=False)
        matched=core.matched_ph_vs_controls(df_all[df_all.geometry.isin(['PH','LINEAR','POLY2','SMOOTHSTEP'])], match_tol=0.35)
        matched.to_csv(resdir/'matched_PH_vs_main_controls_from_stress_tests.csv', index=False)
    else:
        df_all=pd.DataFrame(); matched=pd.DataFrame()
    if env_all:
        pd.concat(env_all, ignore_index=True).to_csv(resdir/'shift_window_envelopes_from_stress_tests.csv', index=False)

    # Plots
    if not df_pert.empty:
        p=df_pert.copy(); p['log10_CF']=np.log10(np.maximum(p['CF_axial_inter_over_horns'],1e-300))
        save_plot_bar(p.sort_values('CF_axial_inter_over_horns', ascending=False), 'geometry', 'CF_axial_inter_over_horns', figdir/'test1_CF_perturbation_controls.png', 'Inter-focal confinement: PHB vs perturbed and control horns', 'CF = E(|x|<=c)/E(|x|>c)', logy=True)
        save_plot_bar(p.sort_values('horn_leakage'), 'geometry', 'horn_leakage', figdir/'test1_horn_leakage_perturbation_controls.png', 'Horn leakage: smaller is better', 'horn leakage fraction', logy=True)
    if not df_cscan.empty:
        save_plot_line(df_cscan, 'c_factor', ['CF_cfactor'], figdir/'test2_CFactor_vs_cfactor.png', 'Where does the axial barrier appear? CF vs c-factor', 'CF(c factor)', logy=True)
        save_plot_line(df_cscan, 'c_factor', ['drop_ratio_cfactor'], figdir/'test2_drop_vs_cfactor.png', 'Focal-plane drop ratio vs diagnostic plane factor', 'drop ratio', logy=True)
    if dens_rows:
        plt.figure(figsize=(9,5))
        for d in dens_rows:
            dd=d
            plt.plot(dd['x'], dd['axial_density'], label=dd['geometry'].iloc[0])
        c=math.sqrt(a*a+b*b)
        plt.axvline(-c, ls=':', color='k'); plt.axvline(c, ls=':', color='k')
        plt.yscale('log'); plt.ylim(bottom=1e-7)
        plt.grid(True, alpha=0.25); plt.legend(); plt.xlabel('x'); plt.ylabel('normalized axial energy density')
        plt.title('Axial density near true focal planes x=±c')
        plt.tight_layout(); plt.savefig(figdir/'test2_axial_density_log.png', dpi=180); plt.close()
    if not df_cap.empty:
        save_plot_bar(df_cap.sort_values('CF_axial_inter_over_horns', ascending=False), 'geometry','CF_axial_inter_over_horns',figdir/'test3_CF_cap_dependence.png','Central cap dependence of inter-focal confinement','CF',logy=True)
    if not df_ms.empty:
        plt.figure(figsize=(9,5))
        for kind,sub in df_ms.groupby('geometry'):
            sub=sub.sort_values('m')
            plt.plot(sub['m'], np.log10(np.maximum(sub['CF_axial_inter_over_horns'],1e-300)), marker='o', label=kind)
        plt.xlabel('m'); plt.ylabel('log10(CF)'); plt.title('m-sweep: axial confinement signature')
        plt.grid(True, alpha=0.25); plt.legend(); plt.tight_layout(); plt.savefig(figdir/'test4_msweep_logCF.png', dpi=180); plt.close()
        plt.figure(figsize=(9,5))
        for kind,sub in df_ms.groupby('geometry'):
            sub=sub.sort_values('m')
            plt.plot(sub['m'], 100*sub['eta_100R'], marker='o', label=kind)
        plt.axhline(70, ls=':', color='k')
        plt.xlabel('m'); plt.ylabel('eta ±10%R, %'); plt.title('m-sweep: annular energy')
        plt.grid(True, alpha=0.25); plt.legend(); plt.tight_layout(); plt.savefig(figdir/'test4_msweep_eta100.png', dpi=180); plt.close()
    if not df_scale.empty:
        plt.figure(figsize=(7,4.5))
        plt.plot(df_scale['scale_factor'], 100*df_scale['eta_100R'], marker='o', label='eta ±10%R')
        plt.plot(df_scale['scale_factor'], np.log10(np.maximum(df_scale['CF_axial_inter_over_horns'],1e-300)), marker='s', label='log10 CF')
        plt.xlabel('scale factor'); plt.title('Scale-invariance smoke test')
        plt.grid(True, alpha=0.25); plt.legend(); plt.tight_layout(); plt.savefig(figdir/'test5_scale_invariance.png', dpi=180); plt.close()
    if not df_confirm.empty:
        save_plot_bar(df_confirm.sort_values('CF_axial_inter_over_horns', ascending=False),'geometry','CF_axial_inter_over_horns',figdir/'test6_confirm100_CF.png','100x70 selected confirmation: CF','CF',logy=True)

    # Text summary.
    summary=[]
    summary.append('PHB C2-C3 FULL UNIQUENESS STRESS TESTS')
    summary.append(f'preset={args.preset}, target a={a}, b={b}, R={R}, m={m_main}, grid main={Nx}x{Ns}')
    summary.append(f'Runtime seconds: {time.time()-t0:.1f}')
    summary.append('')
    if not df_pert.empty:
        ph=df_pert[df_pert.geometry=='PH'].iloc[0]
        summary.append('TEST 1: hyperbolic-law perturbation and controls')
        summary.append(f"PH: eta100={100*ph['eta_100R']:.3f}%, CF={ph['CF_axial_inter_over_horns']:.4e}, horn_leakage={ph['horn_leakage']:.4e}")
        for ctrl in ['LINEAR','POLY2','SMOOTHSTEP','PH_PERT_PLUS_005','PH_PERT_MINUS_005']:
            if not df_pert[df_pert.geometry==ctrl].empty:
                r=df_pert[df_pert.geometry==ctrl].iloc[0]
                summary.append(f"{ctrl}: eta100={100*r['eta_100R']:.3f}%, CF={r['CF_axial_inter_over_horns']:.4e}, PH/ctrl CF={ph['CF_axial_inter_over_horns']/max(r['CF_axial_inter_over_horns'],1e-300):.3g}, ctrl/PH leakage={r['horn_leakage']/max(ph['horn_leakage'],1e-300):.3g}")
        summary.append('')
    if not df_cscan.empty:
        summary.append('TEST 2: focal-plane scan')
        for kind in ['PH','LINEAR','POLY2','SMOOTHSTEP']:
            sub=df_cscan[df_cscan.geometry==kind]
            if not sub.empty:
                mx=sub.loc[sub['CF_cfactor'].idxmax()]
                at1=sub.iloc[(sub['c_factor']-1.0).abs().argsort()[:1]].iloc[0]
                summary.append(f"{kind}: max CF at c_factor={mx['c_factor']:.2f}, max_CF={mx['CF_cfactor']:.4e}; CF at true c={at1['CF_cfactor']:.4e}")
        summary.append('')
    if not df_cap.empty:
        summary.append('TEST 3: cap dependence')
        for _,r in df_cap.sort_values('CF_axial_inter_over_horns',ascending=False).iterrows():
            summary.append(f"{r['geometry']}: eta100={100*r['eta_100R']:.2f}%, CF={r['CF_axial_inter_over_horns']:.3e}, leakage={r['horn_leakage']:.3e}")
        summary.append('')
    if not df_ms.empty:
        summary.append('TEST 4: m-sweep')
        phm=df_ms[df_ms.geometry=='PH'].sort_values('m')
        for _,r in phm.iterrows():
            summary.append(f"PH m={int(r['m'])}: eta100={100*r['eta_100R']:.2f}%, CF={r['CF_axial_inter_over_horns']:.3e}, drop={r['focal_plane_drop_ratio']:.3g}")
        summary.append('')
    if not df_scale.empty:
        summary.append('TEST 5: scale invariance smoke')
        for _,r in df_scale.sort_values('scale_factor').iterrows():
            summary.append(f"scale={r['scale_factor']}: kR={r['kR']:.4f}, eta100={100*r['eta_100R']:.3f}%, CF={r['CF_axial_inter_over_horns']:.3e}")
        summary.append('')
    if not df_confirm.empty:
        summary.append('TEST 6: 100x70 selected confirmation')
        phc=df_confirm[df_confirm.geometry=='PH']
        if not phc.empty:
            ph=phc.iloc[0]
            for _,r in df_confirm.iterrows():
                ratio=ph['CF_axial_inter_over_horns']/max(r['CF_axial_inter_over_horns'],1e-300)
                leakratio=r['horn_leakage']/max(ph['horn_leakage'],1e-300)
                summary.append(f"{r['geometry']}: eta100={100*r['eta_100R']:.3f}%, CF={r['CF_axial_inter_over_horns']:.4e}, PH/this_CF={ratio:.3g}, this/PH_leakage={leakratio:.3g}")
    (resdir/'RUN_SUMMARY_full_uniqueness_stress_tests.txt').write_text('\n'.join(summary), encoding='utf-8')

    # README
    readme=f"""# PHB C2-C3 full uniqueness stress tests

This package contains a reduced scalar Helmholtz stress-test suite for the
closed second-order vertical pseudohyperboloid.  It is aimed at distinguishing
ordinary high-m annular localization from a wave-level imprint of the
hyperbolic focal billiard.

Model status: boundary-fitted axisymmetric scalar Helmholtz eigenproblem.  It
is not a full-vector Maxwell, not a true TE/TM vector calculation, and not a
laser model.

Main script: `PHB_C2C3_full_uniqueness_stress_tests.py`.

Key CSV files are in `results/`; figures are in `figures/`.
"""
    (outdir/'README_REPRODUCIBILITY.md').write_text(readme, encoding='utf-8')

    # Bundle ZIP.
    zip_path=Path(str(outdir)+'_package.zip')
    with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
        z.write(Path(__file__), arcname='PHB_C2C3_full_uniqueness_stress_tests.py')
        z.write(BASE_SCRIPT, arcname='PHB_C2C3_hyperbolic_uniqueness_search_base.py')
        for p in outdir.rglob('*'):
            if p.is_file():
                z.write(p, arcname=str(p.relative_to(outdir)))
    print('\n'.join(summary))
    print('Wrote:', outdir)
    print('ZIP:', zip_path)


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--preset', choices=['quick','standard','confirm'], default='standard')
    ap.add_argument('--outdir', default=str(Path(__file__).resolve().parents[1] / 'run_stress_tests'))
    args=ap.parse_args()
    run_all(args)
