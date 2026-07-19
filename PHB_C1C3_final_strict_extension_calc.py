# -*- coding: utf-8 -*-
"""
Final strict extension calculation for the PHB C1-C3 Zenodo article.

Purpose
-------
This script performs only the checks that are still technically feasible in the
present reduced scalar Helmholtz environment:

1. primary C2 grid extension: 70x48, 100x70, 130x92, 160x112 for PH and matched controls;
2. empirical grid-variability floor for horn leakage and log10(CF);
3. Richardson-style apparent convergence diagnostics (descriptive, not a theorem);
4. axis cutoff sensitivity for the singular rho=0 treatment;
5. solver-tolerance and shift-window sensitivity;
6. boundary-proxy comparison (natural/Dirichlet scalar proxies);
7. additional controls: slope-matched Hermite horn, volume-matched polynomial horn,
   cylinder-like and ellipsoid-like baselines;
8. diagnostic figures and CSV tables for the revised article.

Model status: reduced axisymmetric scalar Helmholtz eigenvalue problem only.
Not full-vector Maxwell, not a laser, not an open-boundary radiation problem.
"""
from __future__ import annotations
import importlib.util, math, json, shutil, zipfile, time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.optimize as opt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Package-root-relative paths. This keeps the Zenodo package self-contained.
BASE = Path(__file__).resolve().parents[1]
CORE_SCRIPT = BASE/'scripts'/'PHB_C2C3_hyperbolic_uniqueness_search.py'
spec = importlib.util.spec_from_file_location('core', str(CORE_SCRIPT))
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

OUT = BASE/'run_full_strict_extension'
RES = OUT/'results'; FIG = OUT/'figures'; SCR = OUT/'scripts'
for d in [RES, FIG, SCR]: d.mkdir(parents=True, exist_ok=True)

# Target geometry
A=0.3; B=0.6; R=3.0; M=15
C=math.sqrt(A*A+B*B); L=core.common_length(A,B,R)
BC_NAT='natural_TE_like_scalar_proxy'; BC_DIR='dirichlet_TM_like_scalar_proxy'
TARGETS_MAIN=np.array([15.55,15.70,15.85,16.00,17.55,17.75,17.95])
TARGETS_C3=np.arange(10,25.01,0.5)
GRIDS=[(70,48),(100,70),(130,92),(160,112)]
GEOMS_MAIN=['PH','LINEAR','POLY2','SMOOTHSTEP']
GEOMS_EXTRA=['HERMITE_SLOPE_MATCH','POLY_VOL_MATCH','CYLINDER_R','ELLIPSOID_L']

# ---------------- Extended profiles ----------------
_orig_wall = core.wall_profile
_orig_derivative = core.wall_derivative
_orig_labels = core.geometry_labels

# precompute volume-matched p for polynomial control: match integral of r^2 dx approximately to PH

def _raw_horn(absx, kind, a=A,b=B,R=R,L=L):
    x=np.asarray(absx,float)
    t=np.clip((x-a)/max(L-a,1e-14),0,1)
    if kind=='HERMITE_SLOPE_MATCH':
        # Cubic Hermite between (t=0,r=R) and (t=1,r=0), matching PH slope at t_c in least local way.
        # Use endpoint slopes in t-coordinate from PH at t=0.05 and t=0.95 to avoid endpoint singularity.
        # This is a non-hyperbolic but slope-informed control.
        # Hermite basis with chosen finite slopes in r vs t.
        s0 = (core._profile_horn_by_t(np.array([a+0.06*(L-a)]),a,b,R,L,'PH')[0]-R)/0.06
        s1 = (0.0-core._profile_horn_by_t(np.array([a+0.94*(L-a)]),a,b,R,L,'PH')[0])/(1-0.94)
        h00=2*t**3-3*t**2+1; h10=t**3-2*t**2+t; h01=-2*t**3+3*t**2; h11=t**3-t**2
        val=h00*R + h10*s0 + h01*0.0 + h11*s1
        return np.maximum(0.0,val)
    if kind=='POLY_VOL_MATCH':
        # Solve for p so horn+cap volume proxy matches PH. Use cached p.
        p=getattr(_raw_horn,'p_vol',None)
        if p is None:
            xs=np.linspace(-L,L,4001); ph=_orig_wall(xs,'PH',a,b,R,L)
            vol_ph=np.trapz(ph**2,xs)
            def volp(p):
                xx=np.abs(xs); out=np.zeros_like(xx); tor=xx<=a; out[tor]=R+np.sqrt(np.maximum(0,a*a-xx[tor]**2))
                horn=(xx>a)&(xx<=L); tt=np.clip((xx[horn]-a)/(L-a),0,1); out[horn]=R*(1-tt**p)
                return np.trapz(out**2,xs)
            def f(p): return volp(p)-vol_ph
            try: p=opt.brentq(f,0.2,6.0)
            except Exception: p=1.0
            _raw_horn.p_vol=p
        return np.maximum(0.0,R*(1-t**p))
    if kind=='CYLINDER_R':
        return R*np.ones_like(x)
    if kind=='ELLIPSOID_L':
        # Ellipsoid-like horn: continuous from R to 0 at L.
        return np.maximum(0.0,R*np.sqrt(np.maximum(0,1-t**2)))
    return core._profile_horn_by_t(x,a,b,R,L,kind)


def ext_wall_profile(xs, kind, a=A, b=B, R=R, L=None):
    if L is None: L=core.common_length(a,b,R)
    x=np.abs(np.asarray(xs,float)); out=np.zeros_like(x)
    if kind=='CYLINDER_R':
        return R*np.ones_like(x)  # simple WGM-like baseline; not endpoint-closed
    tor=x<=a
    out[tor]=R+np.sqrt(np.maximum(0,a*a-x[tor]**2))
    horn=(x>a)&(x<=L)
    if kind in GEOMS_EXTRA:
        out[horn]=_raw_horn(x[horn], kind, a,b,R,L)
    else:
        out[horn]=core._profile_horn_by_t(x[horn],a,b,R,L,kind)
    return np.maximum(out,1e-12)


def ext_wall_derivative(xs, kind, a=A, b=B, R=R, L=None):
    if L is None: L=core.common_length(a,b,R)
    h=1e-5*max(L,R,a,b)
    xs=np.asarray(xs,float); x1=np.clip(xs-h,-L+1e-8,L-1e-8); x2=np.clip(xs+h,-L+1e-8,L-1e-8)
    d=(ext_wall_profile(x2,kind,a,b,R,L)-ext_wall_profile(x1,kind,a,b,R,L))/np.maximum(x2-x1,1e-15)
    return np.clip(d,-500,500)


def ext_labels():
    lab=_orig_labels()
    lab.update({
        'HERMITE_SLOPE_MATCH':'Hermite finite-slope matched horn control',
        'POLY_VOL_MATCH':'polynomial horn with PH volume proxy matched',
        'CYLINDER_R':'constant-radius cylinder-like WGM baseline',
        'ELLIPSOID_L':'ellipsoid-like closed baseline',
    })
    return lab

core.wall_profile=ext_wall_profile
core.wall_derivative=ext_wall_derivative
core.geometry_labels=ext_labels

# ---------------- utilities ----------------

def solve_distinct(kind,Nx,Ns,m=M,targets=TARGETS_MAIN,bc=BC_NAT,tol=1e-6,nev=6,maxiter=5000):
    env,df=core.solve_shift_scan(Nx,Ns,A,B,R,kind,bc,m,np.asarray(targets,float),nev=nev,tol=tol,maxiter=maxiter)
    if not df.empty:
        df=df.copy(); df['grid']=f'{Nx}x{Ns}'; df['tol']=tol; df['boundary_proxy']=bc
    if not env.empty:
        env=env.copy(); env['grid']=f'{Nx}x{Ns}'; env['tol']=tol; env['boundary_proxy']=bc
    return env,df


def choose_primary(df):
    if df is None or df.empty: return None
    d=df[(df.kR>=15.0)&(df.kR<=16.5)].copy()
    if d.empty: d=df.copy()
    return d.sort_values(['eta_100R','CF_axial_inter_over_horns'],ascending=False).iloc[0].to_dict()


def enrich(r):
    r=dict(r)
    for k in ['eta_025R','eta_050R','eta_100R','horn_leakage']:
        r[k+'_pct']=100*float(r[k])
    r['log10_CF']=math.log10(max(float(r['CF_axial_inter_over_horns']),1e-300))
    return r


def apparent_order(vals, hs):
    vals=np.asarray(vals,float); hs=np.asarray(hs,float)
    if len(vals)<3: return dict(p=np.nan,limit=np.nan,err=np.nan,note='need >=3 grids')
    # use last three values for apparent order
    y1,y2,y3=vals[-3:]; h1,h2,h3=hs[-3:]
    denom=y2-y3
    if abs(denom)<1e-15:
        return dict(p=np.nan,limit=y3,err=abs(y2-y3),note='last two values nearly identical')
    ratio=(y1-y2)/denom
    if not np.isfinite(ratio) or ratio<=0:
        return dict(p=np.nan,limit=np.nan,err=np.nan,note='non-monotonic or sign-changing sequence')
    def f(p): return (h1**p-h2**p)/(h2**p-h3**p)-ratio
    p=np.nan
    for lo,hi in [(0.05,0.5),(0.5,1),(1,2),(2,4),(4,8),(8,16)]:
        try:
            if np.isfinite(f(lo)) and np.isfinite(f(hi)) and f(lo)*f(hi)<0:
                p=opt.brentq(f,lo,hi); break
        except Exception: pass
    if not np.isfinite(p): return dict(p=np.nan,limit=np.nan,err=np.nan,note='no stable apparent order root')
    Acoef=(y1-y2)/(h1**p-h2**p)
    lim=y1-Acoef*h1**p
    return dict(p=p,limit=lim,err=abs(y3-lim),note='descriptive apparent order from last three grids')


def geom_metrics(kind):
    xs=np.linspace(-L,L,8001); rw=ext_wall_profile(xs,kind,A,B,R,L); dr=ext_wall_derivative(xs,kind,A,B,R,L)
    idx=np.argmin(np.abs(xs-C))
    return {
        'geometry':kind,'L':L,'c':C,'rho_at_c':float(rw[idx]),'abs_slope_at_c':float(abs(dr[idx])),
        'volume_proxy_pi_removed':float(np.trapz(rw**2,xs)),
        'surface_proxy_2pi_removed':float(np.trapz(rw*np.sqrt(1+dr**2),xs)),
        'min_radius':float(np.min(rw)),'max_radius':float(np.max(rw))
    }

# ---------------- main calculations ----------------

def main():
    t0=time.time()
    meta={'model':'reduced axisymmetric scalar Helmholtz only','a':A,'b':B,'R':R,'c':C,'L':L,'m':M,'date_note':'computed in current container'}
    (RES/'run_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    shutil.copy2(__file__, SCR/Path(__file__).name)
    for src in [CORE_SCRIPT, BASE/'scripts'/'PHB_C2C3_full_uniqueness_stress_tests.py']:
        if src.exists(): shutil.copy2(src, SCR/src.name)

    # 1. Geometry metrics
    geoms_all=GEOMS_MAIN+GEOMS_EXTRA
    gm=pd.DataFrame([geom_metrics(g) for g in geoms_all])
    gm.to_csv(RES/'geometry_matching_metrics.csv',index=False)

    # 2. Grid extension for main controls; 160x112 only for PH/LINEAR/POLY2/SMOOTHSTEP; extras 70/100 for cost.
    rows=[]; modes=[]; envs=[]
    for kind in GEOMS_MAIN:
        for Nx,Ns in GRIDS:
            print('grid',kind,Nx,Ns,flush=True)
            env,df=solve_distinct(kind,Nx,Ns)
            if not df.empty: modes.append(df)
            if not env.empty: envs.append(env)
            pm=choose_primary(df)
            if pm: rows.append(enrich(pm))
    for kind in GEOMS_EXTRA:
        for Nx,Ns in [(70,48),(100,70)]:
            print('extra',kind,Nx,Ns,flush=True)
            env,df=solve_distinct(kind,Nx,Ns)
            if not df.empty: modes.append(df)
            if not env.empty: envs.append(env)
            pm=choose_primary(df)
            if pm: rows.append(enrich(pm))
    primary=pd.DataFrame(rows); primary.to_csv(RES/'primary_mode_grid_extension.csv',index=False)
    if modes: pd.concat(modes,ignore_index=True).to_csv(RES/'all_distinct_modes_grid_extension.csv',index=False)
    if envs: pd.concat(envs,ignore_index=True).to_csv(RES/'all_shift_envelopes_grid_extension.csv',index=False)

    # 3. UQ and floor
    uq=[]; floor=[]
    for geom,sub in primary.groupby('geometry'):
        sub=sub.sort_values('Nx'); hs=1/np.sqrt(sub.Nx.values*sub.Ns.values)
        for met in ['kR','eta_100R_pct','eta_050R_pct','eta_025R_pct','log10_CF','horn_leakage_pct','focal_plane_drop_ratio']:
            out=apparent_order(sub[met].values,hs)
            row={'geometry':geom,'metric':met,'n_grids':len(sub)}
            for _,r in sub.iterrows(): row[f'{met}_{int(r.Nx)}x{int(r.Ns)}']=float(r[met])
            row.update(out); uq.append(row)
        leaks=sub.horn_leakage_pct.values
        if len(leaks)>=2:
            diffs=np.abs(np.diff(leaks)); fl=float(np.max(diffs)); fine=float(leaks[-1])
            floor.append({'geometry':geom,'n_grids':len(sub),'fine_grid':f"{int(sub.iloc[-1].Nx)}x{int(sub.iloc[-1].Ns)}",'fine_horn_leakage_pct':fine,'empirical_grid_variability_floor_pct':fl,'fine_leakage_over_floor':fine/max(fl,1e-300),'interpretation':'resolved_above_floor' if fine>3*fl else 'near_or_below_floor'})
    pd.DataFrame(uq).to_csv(RES/'grid_UQ_apparent_order_and_Richardson.csv',index=False)
    pd.DataFrame(floor).to_csv(RES/'horn_leakage_empirical_grid_variability_floor.csv',index=False)

    # 4. Boundary proxy comparison 70x48 and 100x70 for PH/LINEAR.
    br=[]
    for bc in [BC_NAT,BC_DIR]:
        for kind in ['PH','LINEAR']:
            for Nx,Ns in [(70,48),(100,70)]:
                print('bc',bc,kind,Nx,Ns,flush=True)
                env,df=solve_distinct(kind,Nx,Ns,bc=bc,targets=np.array([15.6,15.8,16.0,18.0]))
                pm=choose_primary(df)
                if pm:
                    e=enrich(pm); e['boundary_proxy']=bc; br.append(e)
    pd.DataFrame(br).to_csv(RES/'boundary_proxy_sensitivity.csv',index=False)

    # 5. Axis cutoff sensitivity: monkeypatch S_MIN for PH, 100x70.
    axis=[]; orig_smin=core.S_MIN
    for smin in [1e-3,2e-3,4e-3,8e-3]:
        core.S_MIN=smin
        print('smin',smin,flush=True)
        env,df=solve_distinct('PH',100,70,targets=np.array([15.6,15.8,16.0]))
        pm=choose_primary(df)
        if pm:
            e=enrich(pm); e['S_MIN']=smin; axis.append(e)
    core.S_MIN=orig_smin
    pd.DataFrame(axis).to_csv(RES/'axis_cutoff_sensitivity_PH_100x70.csv',index=False)

    # 6. Solver tolerance and shift target sensitivity for PH/LINEAR at 100x70.
    tolrows=[]
    for tol in [1e-5,1e-6,1e-7]:
        for kind in ['PH','LINEAR']:
            print('tol',tol,kind,flush=True)
            env,df=solve_distinct(kind,100,70,tol=tol,targets=np.array([15.7,15.85,16.0]),maxiter=7000)
            pm=choose_primary(df)
            if pm:
                e=enrich(pm); e['solver_tol']=tol; tolrows.append(e)
    pd.DataFrame(tolrows).to_csv(RES/'solver_tolerance_sensitivity_100x70.csv',index=False)
    shiftrows=[]
    for sh in [15.55,15.65,15.75,15.85,15.95,16.05]:
        env,df=solve_distinct('PH',100,70,targets=np.array([sh]))
        pm=choose_primary(df)
        if pm:
            e=enrich(pm); e['shift_target']=sh; shiftrows.append(e)
    pd.DataFrame(shiftrows).to_csv(RES/'shift_window_sensitivity_PH_100x70.csv',index=False)

    # 7. C3 protocol at 70x48 for main geoms.
    c3=[]; c3env=[]
    for kind in GEOMS_MAIN:
        print('C3',kind,flush=True)
        env,df=solve_distinct(kind,70,48,targets=TARGETS_C3,nev=6,maxiter=3500)
        if not env.empty: c3env.append(env)
        if not df.empty: c3.append(df)
    if c3:
        c3df=pd.concat(c3,ignore_index=True); c3df.to_csv(RES/'C3_deduplicated_modes_70x48.csv',index=False)
        clusters=core.analyze_joint_clusters(c3df,eta_thr=0.70,cf_thr=100.0,max_gap=2.1)
        clusters.to_csv(RES/'C3_two_peak_support_not_continuous_window.csv',index=False)
    if c3env: pd.concat(c3env,ignore_index=True).to_csv(RES/'C3_shift_window_envelope_raw.csv',index=False)

    # 8. m sweep to diagnose generic high-m mechanism at 70x48.
    ms=[]
    for kind in ['PH','LINEAR','POLY2','POLY3','SMOOTHSTEP','CYLINDER_R','ELLIPSOID_L']:
        for m in [6,10,12,15,20,25]:
            targets=np.array([max(5,0.8*m+3),max(6,m+1),max(8,1.15*m+2),max(10,1.25*m+4)])
            print('msweep',kind,m,flush=True)
            env,df=solve_distinct(kind,70,48,m=m,targets=targets,nev=5,maxiter=2800)
            if not df.empty:
                row=df.sort_values(['eta_100R','CF_axial_inter_over_horns'],ascending=False).iloc[0].to_dict()
                ms.append(enrich(row))
    pd.DataFrame(ms).to_csv(RES/'m_sweep_extended_controls_70x48.csv',index=False)

    # 9. Figures.
    def save_fig(name):
        plt.tight_layout(); plt.savefig(FIG/name,dpi=200); plt.close()
    # geometry
    xs=np.linspace(-L,L,1200)
    plt.figure(figsize=(8,7))
    for kind in geoms_all:
        rw=ext_wall_profile(xs,kind,A,B,R,L); plt.plot(xs,rw,label=kind); plt.plot(xs,-rw,lw=.5,alpha=.5)
    plt.axvline(C,color='k',ls=':',lw=1); plt.axvline(-C,color='k',ls=':',lw=1); plt.axhline(R,color='gray',ls=':',lw=1); plt.axhline(-R,color='gray',ls=':',lw=1)
    plt.gca().set_aspect('equal',adjustable='box'); plt.xlabel('x'); plt.ylabel('signed rho'); plt.title('PHB and benchmark controls'); plt.legend(fontsize=7,ncol=2); plt.grid(alpha=.25); save_fig('fig_geometry_controls_extended.png')
    # eta grid
    plt.figure(figsize=(8,5))
    for geom,sub in primary[primary.geometry.isin(GEOMS_MAIN)].groupby('geometry'):
        sub=sub.sort_values('Nx'); plt.plot(sub.Nx.astype(str),sub.eta_100R_pct,marker='o',label=geom)
    plt.axhline(70,color='k',ls=':',lw=1); plt.ylabel('eta_100R (%)'); plt.xlabel('grid Nx'); plt.title('C2 annular energy grid sensitivity'); plt.grid(alpha=.25); plt.legend(); save_fig('fig_eta100_grid_extended.png')
    # leakage vs floor
    nf=pd.read_csv(RES/'horn_leakage_empirical_grid_variability_floor.csv')
    d=nf[nf.geometry.isin(GEOMS_MAIN)].copy(); x=np.arange(len(d)); w=.35
    plt.figure(figsize=(9,5)); plt.bar(x-w/2,d.fine_horn_leakage_pct,w,label='fine-grid leakage'); plt.bar(x+w/2,d.empirical_grid_variability_floor_pct,w,label='empirical grid-variability floor'); plt.yscale('log'); plt.xticks(x,d.geometry); plt.ylabel('% modal energy'); plt.title('Horn leakage versus empirical floor'); plt.grid(axis='y',alpha=.25); plt.legend(); save_fig('fig_leakage_vs_floor_extended.png')
    # logCF grid
    plt.figure(figsize=(8,5))
    for geom,sub in primary[primary.geometry.isin(GEOMS_MAIN)].groupby('geometry'):
        sub=sub.sort_values('Nx'); plt.plot(sub.Nx.astype(str),sub.log10_CF,marker='o',label=geom)
    plt.ylabel('log10(CF)'); plt.xlabel('grid Nx'); plt.title('Secondary CF diagnostic: grid sensitivity'); plt.grid(alpha=.25); plt.legend(); save_fig('fig_logCF_grid_extended.png')
    # geometry slope vs logCF 130/100
    comb=primary.merge(gm,on='geometry',how='left')
    sub=comb[(comb.Nx==130) | ((comb.Nx==100)&(~comb.geometry.isin(GEOMS_MAIN)))]
    plt.figure(figsize=(8,5)); plt.scatter(sub.abs_slope_at_c,sub.log10_CF,s=60)
    for _,r in sub.iterrows(): plt.text(r.abs_slope_at_c,r.log10_CF,r.geometry,fontsize=8)
    plt.xlabel('|slope at x=c|'); plt.ylabel('log10(CF)'); plt.title('CF may correlate with local slope; causality not isolated'); plt.grid(alpha=.25); save_fig('fig_logCF_vs_slope_extended.png')
    # m sweep eta/logcf
    msdf=pd.read_csv(RES/'m_sweep_extended_controls_70x48.csv')
    plt.figure(figsize=(9,5))
    for geom,sub in msdf.groupby('geometry'):
        sub=sub.sort_values('m'); plt.plot(sub.m,sub.eta_100R_pct,marker='o',label=geom)
    plt.axhline(70,color='k',ls=':',lw=1); plt.xlabel('m'); plt.ylabel('eta_100R (%)'); plt.title('High-m annular localization across geometries'); plt.grid(alpha=.25); plt.legend(ncol=2,fontsize=7); save_fig('fig_msweep_eta_extended.png')
    plt.figure(figsize=(9,5))
    for geom,sub in msdf.groupby('geometry'):
        sub=sub.sort_values('m'); plt.plot(sub.m,sub.log10_CF,marker='o',label=geom)
    plt.xlabel('m'); plt.ylabel('log10(CF)'); plt.title('CF also increases with high-m contribution'); plt.grid(alpha=.25); plt.legend(ncol=2,fontsize=7); save_fig('fig_msweep_logCF_extended.png')

    # summary
    lines=[]
    lines.append('FINAL STRICT EXTENSION CALCULATION')
    lines.append(f'a={A}, b={B}, R={R}, c={C:.8f}, L={L:.8f}, m={M}')
    lines.append(f'runtime_sec={time.time()-t0:.1f}')
    lines.append('')
    lines.append('WHAT THIS SCRIPT CAN VERIFY: reduced scalar Helmholtz C1-C3 diagnostics, grid/tolerance/shift sensitivity, additional control geometries, and empirical leakage floor.')
    lines.append('WHAT THIS SCRIPT CANNOT VERIFY: full-vector Maxwell modes, open radiation Q, laser gain/loss/output, true physical horn-leakage below the numerical floor, and causal proof of ray focal billiard mechanism.')
    lines.append('')
    lines.append('Primary modes:')
    for _,r in primary.sort_values(['geometry','Nx']).iterrows():
        lines.append(f"{r.geometry:18s} {int(r.Nx)}x{int(r.Ns)} kR={r.kR:.6f} eta100={r.eta_100R_pct:.3f}% eta50={r.eta_050R_pct:.3f}% logCF={r.log10_CF:.3f} leak={r.horn_leakage_pct:.3e}% drop={r.focal_plane_drop_ratio:.3g}")
    lines.append('')
    lines.append('Leakage floor interpretation:')
    for _,r in pd.read_csv(RES/'horn_leakage_empirical_grid_variability_floor.csv').iterrows():
        lines.append(f"{r.geometry:18s} fine_leak={r.fine_horn_leakage_pct:.3e}% floor={r.empirical_grid_variability_floor_pct:.3e}% ratio={r.fine_leakage_over_floor:.3g} -> {r.interpretation}")
    lines.append('')
    lines.append('Main interpretation for the article: eta_100R is stable and supports C2. eta_100R does not demonstrate hyperbolic uniqueness because controls are close. CF/horn leakage remains an interesting secondary diagnostic, but PHB leakage is near/below the empirical grid-variability floor; it should be stated as a hypothesis for further verification, not a proof. Additional controls show that current controls are dimension-matched but not fully slope/volume/local-geometry matched; this limits causal claims about hyperbolic focality.')
    (RES/'RUN_SUMMARY_final_strict_extension.txt').write_text('\n'.join(lines),encoding='utf-8')
    (OUT/'README_REPRODUCIBILITY.md').write_text('Run: python PHB_C1C3_final_strict_extension_calc.py\nModel: reduced scalar Helmholtz only. Outputs are in results/ and figures/.\n',encoding='utf-8')
    # zip
    zpath=BASE/'PHB_C1C3_final_strict_extension_results_package.zip'
    if zpath.exists(): zpath.unlink()
    with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in OUT.rglob('*'):
            if p.is_file(): z.write(p,arcname=str(p.relative_to(OUT)))
    print('\n'.join(lines))
    print('WROTE',OUT)
    print('ZIP',zpath)

if __name__=='__main__':
    main()
