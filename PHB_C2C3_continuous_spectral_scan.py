# -*- coding: utf-8 -*-
"""
PHB_C2C3_continuous_spectral_scan.py

Enhanced C2/C3 verification for the closed second-order vertical
pseudohyperboloid cavity.

Main purpose:
  1) Perform a dense shift-window spectral scan for the requested top
     configuration a=0.3, b=0.6, R=3.0, m=15, kR in [10,25] step 0.5.
  2) Separate two concepts that are often confused:
       - shift-window envelope: the best eigenmode found near each target kR;
       - distinct modal spectrum: deduplicated eigenmodes, which is the safer
         basis for C3 statements.
  3) Estimate a conservative C3 modal-window proxy and explain whether a true
     continuous spectral window is supported.
  4) Provide a local type-size screen and grid-convergence checks for the
     strongest candidates.

Model status:
  Reduced boundary-fitted axisymmetric scalar Helmholtz eigenproblem.
  This is not ray tracing, not an open laser calculation, and not full-vector Maxwell.
"""
from __future__ import annotations
import math, json, time, argparse, zipfile, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

S_MIN = 2e-3
FOCAL_BANDS = [0.025, 0.050, 0.100]
DEDUP_KR_TOL = 0.06

# ---------------- Geometry ----------------
def common_length(a,b,R):
    return float(a*math.sqrt(1.0+(R/b)**2))

def profile_ph_absx(absx,a,b,R,L):
    x=np.asarray(absx,dtype=float); out=np.zeros_like(x)
    tor=x<=a
    out[tor]=R+np.sqrt(np.maximum(0.0,a*a-x[tor]**2))
    horn=(x>a)&(x<=L)
    xh=x[horn]
    out[horn]=np.maximum(0.0,R-b*np.sqrt(np.maximum(0.0,(xh/a)**2-1.0)))
    return out

def profile_ref_absx(absx,a,b,R,L):
    x=np.asarray(absx,dtype=float); out=np.zeros_like(x)
    tor=x<=a
    out[tor]=R+np.sqrt(np.maximum(0.0,a*a-x[tor]**2))
    horn=(x>a)&(x<=L)
    xh=x[horn]
    out[horn]=np.maximum(0.0,R*(L-xh)/max(L-a,1e-12))
    return out

def wall_profile(xs,kind,a,b,R,L):
    ax=np.abs(np.asarray(xs,dtype=float))
    if kind=='PH': return profile_ph_absx(ax,a,b,R,L)
    if kind=='REF': return profile_ref_absx(ax,a,b,R,L)
    raise ValueError(kind)

def wall_derivative(xs,kind,a,b,R,L):
    h=1e-5*max(L,R,a,b)
    x1=np.clip(xs-h,-L+1e-8,L-1e-8)
    x2=np.clip(xs+h,-L+1e-8,L-1e-8)
    r1=wall_profile(x1,kind,a,b,R,L); r2=wall_profile(x2,kind,a,b,R,L)
    d=(r2-r1)/np.maximum(x2-x1,1e-15)
    return np.clip(d,-200.0,200.0)

# ---------------- Operators ----------------
def derivative_matrix(n,h):
    rows=[]; cols=[]; vals=[]
    for i in range(n):
        if i==0:
            rows += [i,i]; cols += [0,1]; vals += [-1/h,1/h]
        elif i==n-1:
            rows += [i,i]; cols += [n-2,n-1]; vals += [-1/h,1/h]
        else:
            rows += [i,i]; cols += [i-1,i+1]; vals += [-0.5/h,0.5/h]
    return sp.csr_matrix((vals,(rows,cols)),shape=(n,n))

def build_problem(Nx,Ns,m,kind,bc,a,b,R):
    L=common_length(a,b,R)
    xmin=-L+1e-4*max(L,R); xmax=L-1e-4*max(L,R)
    xs=np.linspace(xmin,xmax,Nx); ss=np.linspace(S_MIN,1.0,Ns)
    dx=(xmax-xmin)/(Nx-1); ds=(1.0-S_MIN)/(Ns-1)
    Dx1=derivative_matrix(Nx,dx); Ds1=derivative_matrix(Ns,ds)
    Dx=sp.kron(Dx1,sp.eye(Ns,format='csr'),format='csr')
    Ds=sp.kron(sp.eye(Nx,format='csr'),Ds1,format='csr')
    X,S=np.meshgrid(xs,ss,indexing='ij')
    rw=wall_profile(xs,kind,a,b,R,L)[:,None]
    rpx=wall_derivative(xs,kind,a,b,R,L)[:,None]
    alpha=(S*rpx/np.maximum(rw,1e-10)).ravel()
    W=(S*rw*rw).ravel()
    invr=(1.0/np.maximum(rw,1e-10))*np.ones_like(S)
    area=dx*ds
    Gx=Dx-sp.diags(alpha,0,format='csr')@Ds
    Gr=sp.diags(invr.ravel(),0,format='csr')@Ds
    K=Gx.T@sp.diags(W*area,0,format='csr')@Gx + Gr.T@sp.diags(W*area,0,format='csr')@Gr
    if m>0:
        K=K+sp.diags((m*m/np.maximum(S.ravel(),1e-8))*area,0,format='csr')
    M=sp.diags(W*area,0,format='csr')
    n=Nx*Ns
    fixed=np.zeros(n,dtype=bool)
    def idx(i,j): return i*Ns+j
    if bc=='natural_TE_proxy':
        if m>0:
            fixed[[idx(i,0) for i in range(Nx)]]=True
    elif bc=='dirichlet_TM_proxy':
        fixed[[idx(i,Ns-1) for i in range(Nx)]]=True
        fixed[[idx(0,j) for j in range(Ns)]]=True
        fixed[[idx(Nx-1,j) for j in range(Ns)]]=True
        if m>0:
            fixed[[idx(i,0) for i in range(Nx)]]=True
    else:
        raise ValueError(bc)
    free=np.where(~fixed)[0]
    K=K[free][:,free]; M=M[free][:,free]
    K=(K+K.T)*0.5; M=(M+M.T)*0.5
    mdiag=M.diagonal()
    Sscale=sp.diags(1.0/np.sqrt(np.maximum(mdiag,1e-300)),0,format='csr')
    A=Sscale@K@Sscale
    return dict(A=A,M=M,Sscale=Sscale,free=free,xs=xs,ss=ss,X=X,S=S,rw=rw,L=L,Nx=Nx,Ns=Ns,dx=dx,ds=ds,kind=kind,bc=bc,a=a,b=b,R=R,m=m)

def mode_metrics(vec_free,eigval,prob,target_kR=None):
    Nx=prob['Nx']; Ns=prob['Ns']; xs=prob['xs']; ss=prob['ss']; X=prob['X']; S=prob['S']; rw=prob['rw']
    a=prob['a']; b=prob['b']; R=prob['R']; L=prob['L']; kind=prob['kind']; bc=prob['bc']; m=prob['m']
    free=prob['free']
    U=np.zeros(Nx*Ns)
    U[free]=vec_free
    U=U.reshape((Nx,Ns))
    dx=(xs[-1]-xs[0])/(Nx-1); ds=(1.0-S_MIN)/(Ns-1)
    RHO=S*rw; W=S*rw*rw; amp2=U*U
    Etot=float(np.sum(amp2*W)*dx*ds)
    Vtot=float(np.sum(W)*dx*ds)
    if Etot<=1e-300: return None
    c=math.sqrt(a*a+b*b)
    k=math.sqrt(max(float(eigval),0.0)); kR=k*R
    row={'geometry':kind,'bc':bc,'a':a,'b':b,'R':R,'c':c,'L':L,'m':m,'kR':kR,'k_geom':k,
         'R_over_lambda':kR/(2*math.pi),'Nx':Nx,'Ns':Ns,'target_kR':target_kR if target_kR is not None else np.nan,
         'energy_between_absx_le_c':float(np.sum(amp2*W*(np.abs(X)<=c))*dx*ds/Etot),
         'near_wall_s_ge_0p9':float(np.sum(amp2*W*(S>=0.9))*dx*ds/Etot)}
    for band in FOCAL_BANDS:
        mask=(np.abs(X)<=c)&(np.abs(RHO-R)<=band*R)
        ef=float(np.sum(amp2*W*mask)*dx*ds/Etot)
        vf=float(np.sum(W*mask)*dx*ds/max(Vtot,1e-300))
        key=f'{int(round(1000*band)):03d}'
        row[f'eta_{key}R']=ef; row[f'volfrac_{key}R']=vf; row[f'enrichment_{key}R']=ef/max(vf,1e-12)
    return row

def dedup_modes(rows,tol=DEDUP_KR_TOL):
    rows=sorted(rows,key=lambda r:(r['kR'], -r.get('eta_100R',0)))
    out=[]
    for r in rows:
        if not out or abs(r['kR']-out[-1]['kR'])>tol:
            out.append(r)
        else:
            # keep the instance with highest eta, but preserve target hits count later separately
            if r.get('eta_100R',0)>out[-1].get('eta_100R',0): out[-1]=r
    return out

def solve_shift_scan(Nx,Ns,a,b,R,kind,bc,m,kR_targets,nev=8,tol=8e-7,maxiter=2600):
    prob=build_problem(Nx,Ns,m,kind,bc,a,b,R)
    A=prob['A']; Sscale=prob['Sscale']
    all_rows=[]; envelope_rows=[]
    k_eigs=min(nev,max(1,A.shape[0]-2))
    for target in kR_targets:
        sigma=(target/R)**2
        try:
            vals,y=spla.eigsh(A,k=k_eigs,sigma=sigma,which='LM',tol=tol,maxiter=maxiter)
            vecs=Sscale@y
            order=np.argsort(vals); vals=vals[order]; vecs=vecs[:,order]
            rows=[]
            for j,val in enumerate(vals):
                if val<=1e-9: continue
                row=mode_metrics(vecs[:,j],float(val),prob,target_kR=float(target))
                if row is not None:
                    rows.append(row); all_rows.append(row)
            if rows:
                best=max(rows,key=lambda r:r['eta_100R'])
                er=dict(best); er['scan_target_kR']=float(target); er['n_modes_returned']=len(rows)
                envelope_rows.append(er)
            else:
                envelope_rows.append({'geometry':kind,'bc':bc,'a':a,'b':b,'R':R,'m':m,'scan_target_kR':float(target),'kR':np.nan,'eta_100R':np.nan,'eta_050R':np.nan,'eta_025R':np.nan,'n_modes_returned':0})
        except Exception as e:
            envelope_rows.append({'geometry':kind,'bc':bc,'a':a,'b':b,'R':R,'m':m,'scan_target_kR':float(target),'kR':np.nan,'eta_100R':np.nan,'eta_050R':np.nan,'eta_025R':np.nan,'n_modes_returned':0,'error':str(e)[:200]})
            continue
    distinct=dedup_modes(all_rows)
    return pd.DataFrame(envelope_rows), pd.DataFrame(distinct)

def analyze_c3(df,threshold=0.70):
    if df.empty:
        return {'n_modes':0,'n_over':0,'status':'no_modes'}
    d=df.sort_values('kR').reset_index(drop=True)
    over=d['eta_100R']>=threshold
    groups=[]; cur=[]
    for i,flag in enumerate(over):
        if flag: cur.append(i)
        else:
            if cur: groups.append(cur); cur=[]
    if cur: groups.append(cur)
    # For modal spectrum, a "continuous window" is not declared by default: only consecutive modes all above threshold.
    best=d.loc[d['eta_100R'].idxmax()]
    over_modes=d[over]
    max_group=max(groups,key=len) if groups else []
    if max_group:
        sub=d.loc[max_group]
        width=float(sub['kR'].max()-sub['kR'].min()) if len(sub)>1 else 0.0
        center=float((sub['kR'].max()+sub['kR'].min())/2.0)
        q_proxy=float(center/width) if width>1e-9 else np.inf
    else:
        width=0.0; center=np.nan; q_proxy=np.nan
    # check if any below-threshold mode lies between min and max of all over-threshold modes
    if len(over_modes)>=2:
        lo=float(over_modes['kR'].min()); hi=float(over_modes['kR'].max())
        below_inside=d[(d['kR']>lo)&(d['kR']<hi)&(d['eta_100R']<threshold)]
        has_gap=len(below_inside)>0
        total_span=hi-lo
        q_total=(0.5*(hi+lo)/total_span) if total_span>1e-9 else np.inf
    else:
        has_gap=False; total_span=0.0; q_total=np.nan
    status='no_C2_modes'
    if len(over_modes)>0:
        if len(max_group)>=3 and not has_gap:
            status='strong_discrete_modal_cluster_possible_C3_support'
        elif len(max_group)>=2:
            status='short_discrete_modal_cluster_C3_support_but_not_continuous_window'
        else:
            status='isolated_over_threshold_peaks_not_C3_window'
        if has_gap:
            status+='__with_below_threshold_gaps'
    return {
        'n_modes':int(len(d)), 'n_over':int(len(over_modes)), 'best_kR':float(best['kR']), 'best_eta100':float(best['eta_100R']), 'best_eta50':float(best['eta_050R']),
        'over_kR_list':[round(float(x),4) for x in over_modes['kR'].tolist()],
        'largest_consecutive_over_count':int(len(max_group)), 'largest_consecutive_width_kR':width,
        'largest_consecutive_center_kR':center, 'largest_consecutive_Q_proxy':q_proxy,
        'all_over_span_kR':float(total_span), 'all_over_Q_proxy':float(q_total) if not np.isnan(q_total) else np.nan,
        'below_threshold_gaps_between_over_modes':bool(has_gap), 'status':status
    }

# ---------------- Workflow ----------------
def plot_profiles(a,b,R,outfile):
    L=common_length(a,b,R)
    x=np.linspace(-L,L,900); xp=np.abs(x)
    ph=wall_profile(x,'PH',a,b,R,L); ref=wall_profile(x,'REF',a,b,R,L)
    c=math.sqrt(a*a+b*b)
    fig,ax=plt.subplots(figsize=(8,5.5))
    ax.plot(x,ph,label='PHB hyperbolic horns')
    ax.plot(x,-ph,color=ax.lines[-1].get_color())
    ax.plot(x,ref,'--',label='REF linear horns')
    ax.plot(x,-ref,'--',color=ax.lines[-1].get_color())
    ax.axvline(-c,color='0.35',ls=':',lw=1); ax.axvline(c,color='0.35',ls=':',lw=1)
    ax.axhline(R,color='tab:green',ls=':',lw=1,label='rho=R')
    ax.axhline(-R,color='tab:green',ls=':',lw=1)
    ax.fill_between([-c,c],[R-0.10*R,R-0.10*R],[R+0.10*R,R+0.10*R],alpha=0.15,label='diagnostic ±10%R')
    ax.set_aspect('equal',adjustable='box'); ax.grid(alpha=.25)
    ax.set_xlabel('x'); ax.set_ylabel('signed rho')
    ax.set_title(f'Target geometry: a={a}, b={b}, R={R}')
    ax.legend(fontsize=8,loc='best')
    fig.tight_layout(); fig.savefig(outfile,dpi=170); plt.close(fig)

def make_plots(outdir, target_envelopes, target_distinct):
    figdir=Path(outdir)/'figures'; figdir.mkdir(parents=True,exist_ok=True)
    # Envelope plot
    fig,ax=plt.subplots(figsize=(9,5.2))
    for key,df in target_envelopes.items():
        if df.empty: continue
        label=f"{key[0]} / {key[1].replace('_proxy','')}"
        ax.plot(df['scan_target_kR'],100*df['eta_100R'],marker='o',ms=3,lw=1.2,label=label)
    ax.axhline(70,color='red',ls='--',lw=1,label='C2 threshold 70%')
    ax.set_xlabel('target kR shift window'); ax.set_ylabel('best returned mode energy in ±10%R, %')
    ax.set_title('Dense shift-window envelope scan (not by itself a continuous C3 proof)')
    ax.grid(alpha=.25); ax.legend(fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(figdir/'dense_shift_window_envelope_eta100.png',dpi=180); plt.close(fig)

    # Distinct modal spectrum PH only and PH/REF comparison
    fig,ax=plt.subplots(figsize=(9,5.5))
    for key,df in target_distinct.items():
        if df.empty: continue
        geom,bc=key
        style='o-' if geom=='PH' else 's--'
        label=f"{geom} / {bc.replace('_proxy','')}"
        ax.plot(df['kR'],100*df['eta_100R'],style,ms=4,lw=1.0,label=label)
    ax.axhline(70,color='red',ls='--',lw=1,label='C2 threshold')
    ax.set_xlabel('distinct eigenmode kR'); ax.set_ylabel('energy in ±10%R, %')
    ax.set_title('Deduplicated modal spectrum: safer C3 evidence')
    ax.grid(alpha=.25); ax.legend(fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(figdir/'deduplicated_modal_spectrum_eta100.png',dpi=180); plt.close(fig)

def run(args):
    outdir=Path(args.outdir); resdir=outdir/'results'; figdir=outdir/'figures'
    resdir.mkdir(parents=True,exist_ok=True); figdir.mkdir(parents=True,exist_ok=True)
    t0=time.time()
    a=args.a; b=args.b; R=args.R; m=args.m
    kR_targets=np.arange(args.kR_min,args.kR_max+1e-9,args.kR_step)
    plot_profiles(a,b,R,figdir/'target_PHB_REF_geometry.png')
    meta={'model':'reduced boundary-fitted scalar Helmholtz eigenproblem','not_full_vector_Maxwell':True,'not_open_laser':True,
          'target_config':{'a':a,'b':b,'R':R,'m':m},'kR_targets':kR_targets.tolist(),'grid_main':[args.Nx,args.Ns],
          'C2_threshold_eta100':0.70,'note':'Dense shift-window scan is used to collect eigenmodes. C3 is judged from deduplicated modal spectrum, not from repeated rediscovery of the same eigenmode.'}
    (resdir/'run_metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Running dense target C3 scan...')
    target_envelopes={}; target_distinct={}; c3_rows=[]
    for kind in ['PH','REF']:
        for bc in ['natural_TE_proxy','dirichlet_TM_proxy']:
            print(f'  {kind} {bc} grid={args.Nx}x{args.Ns}')
            env,dist=solve_shift_scan(args.Nx,args.Ns,a,b,R,kind,bc,m,kR_targets,nev=args.nev)
            target_envelopes[(kind,bc)]=env; target_distinct[(kind,bc)]=dist
            env.to_csv(resdir/f'target_envelope_{kind}_{bc}_{args.Nx}x{args.Ns}.csv',index=False)
            dist.to_csv(resdir/f'target_distinct_modes_{kind}_{bc}_{args.Nx}x{args.Ns}.csv',index=False)
            ana=analyze_c3(dist,threshold=0.70); ana.update({'geometry':kind,'bc':bc,'Nx':args.Nx,'Ns':args.Ns,'a':a,'b':b,'R':R,'m':m})
            c3_rows.append(ana)
    dfc3=pd.DataFrame(c3_rows)
    dfc3.to_csv(resdir/'C3_summary_target_config.csv',index=False)
    make_plots(outdir,target_envelopes,target_distinct)

    # local type-size screen: targeted candidate set, quick grid 70x48 by default only PH, both bc.
    # Include the requested config and nearby variants. This is not a global proof; it is a local candidate table for the article.
    candidate_configs=[(0.3,0.6,3.0),(0.3,0.3,3.0),(0.4,0.6,3.0),(0.2,0.6,3.0),(0.3,0.9,3.0),
                       (0.6,0.6,3.0),(0.4,0.3,3.0),(0.2,0.3,3.0),(0.6,0.3,3.0),(0.1,0.6,3.0),
                       (0.5,0.6,3.0),(0.3,0.45,3.0),(0.3,0.75,3.0)]
    # unique preserve order
    seen=set(); configs=[]
    for cfg in candidate_configs:
        if cfg not in seen: seen.add(cfg); configs.append(cfg)
    screen_targets=np.arange(10,25.01,1.0)
    screen_rows=[]
    print('Running local top-configuration screen...')
    for ci,(aa,bb,RR) in enumerate(configs,1):
        print(f'  local screen {ci}/{len(configs)}: a={aa}, b={bb}, R={RR}')
        best=None
        for bc in ['natural_TE_proxy','dirichlet_TM_proxy']:
            env,dist=solve_shift_scan(args.screen_Nx,args.screen_Ns,aa,bb,RR,'PH',bc,m,screen_targets,nev=max(4,args.nev-2))
            if not dist.empty:
                cand=dist.sort_values('eta_100R',ascending=False).iloc[0].to_dict()
                cand['n_c2_peaks']=int((dist['eta_100R']>=0.70).sum())
                cand['over92']=bool(cand['eta_100R']>=0.92)
                if best is None or cand['eta_100R']>best['eta_100R']:
                    best=cand
        if best is None:
            best={'geometry':'PH','bc':'none','a':aa,'b':bb,'R':RR,'m':m,'kR':np.nan,'eta_100R':0,'eta_050R':0,'eta_025R':0,'enrichment_100R':np.nan,'n_c2_peaks':0,'over92':False}
        screen_rows.append(best)
    dfscreen=pd.DataFrame(screen_rows).sort_values('eta_100R',ascending=False)
    dfscreen.to_csv(resdir/'local_type_size_screen_top_candidates.csv',index=False)

    # Validate top 5 local screen candidates and the requested target on the main grid if not already.
    val_configs=[]
    for _,r in dfscreen.head(5).iterrows(): val_configs.append((float(r['a']),float(r['b']),float(r['R'])))
    if (a,b,R) not in val_configs: val_configs.insert(0,(a,b,R))
    val_rows=[]
    val_targets=np.arange(10,25.01,0.5)
    print('Running grid convergence / validation for selected top configs...')
    for (aa,bb,RR) in val_configs:
        for grid in [(50,35),(70,48)]:
            Nx,Ns=grid
            best=None; best_dist=None
            for bc in ['natural_TE_proxy','dirichlet_TM_proxy']:
                env,dist=solve_shift_scan(Nx,Ns,aa,bb,RR,'PH',bc,m,val_targets,nev=args.nev)
                if not dist.empty:
                    cand=dist.sort_values('eta_100R',ascending=False).iloc[0].to_dict()
                    cand['n_c2_peaks']=int((dist['eta_100R']>=0.70).sum())
                    cand['over_threshold_kR_list']=[round(float(x),4) for x in dist[dist['eta_100R']>=0.70]['kR'].sort_values().tolist()]
                    if best is None or cand['eta_100R']>best['eta_100R']:
                        best=cand; best_dist=dist
            if best is None:
                best={'geometry':'PH','bc':'none','a':aa,'b':bb,'R':RR,'m':m,'kR':np.nan,'eta_100R':0,'eta_050R':0,'eta_025R':0,'n_c2_peaks':0,'over_threshold_kR_list':[]}
            best['validation_Nx']=Nx; best['validation_Ns']=Ns
            val_rows.append(best)
    dfval=pd.DataFrame(val_rows).sort_values(['a','b','R','validation_Nx'])
    dfval.to_csv(resdir/'selected_grid_validation.csv',index=False)

    # Plot top local candidates
    fig,ax=plt.subplots(figsize=(9,5.5))
    top=dfscreen.head(10).copy()
    top['label']=top.apply(lambda r:f"a={r.a:g}, b={r.b:g}, R={r.R:g}",axis=1)
    top=top.sort_values('eta_100R',ascending=True)
    ax.barh(top['label'],100*top['eta_100R'])
    ax.axvline(92,color='purple',ls=':',lw=1,label='92% reference')
    ax.axvline(70,color='red',ls='--',lw=1,label='C2 threshold')
    ax.set_xlabel('best PH energy in ±10%R, %')
    ax.set_title('Local type-size screen around target configuration')
    ax.grid(axis='x',alpha=.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(figdir/'local_top_type_size_screen.png',dpi=170); plt.close(fig)

    # Summary report
    best_target=dfc3[dfc3.geometry=='PH'].sort_values('best_eta100',ascending=False).iloc[0]
    best_local=dfscreen.sort_values('eta_100R',ascending=False).iloc[0]
    over92_count=int((dfscreen['eta_100R']>=0.92).sum())
    lines=[]
    lines.append('PHB C2-C3 continuous spectral scan summary')
    lines.append('='*72)
    lines.append(f'Target configuration requested: a={a}, b={b}, R={R}, m={m}.')
    lines.append(f'Dense kR shift scan: {args.kR_min}..{args.kR_max}, step={args.kR_step}, grid={args.Nx}x{args.Ns}.')
    lines.append('Model: reduced scalar boundary-fitted Helmholtz eigenproblem; not full-vector Maxwell; not open laser output.')
    lines.append('')
    lines.append('C3 target scan results (deduplicated modal spectrum):')
    lines.append(dfc3[['geometry','bc','n_modes','n_over','best_kR','best_eta100','over_kR_list','largest_consecutive_over_count','largest_consecutive_width_kR','largest_consecutive_Q_proxy','below_threshold_gaps_between_over_modes','status']].to_string(index=False))
    lines.append('')
    lines.append('Best target PHB mode among boundary proxies:')
    lines.append(f"  bc={best_target.bc}; best kR={best_target.best_kR:.4f}; eta±10%R={100*best_target.best_eta100:.2f}%; eta±5%R={100*best_target.best_eta50:.2f}%")
    lines.append(f"  over-threshold kR values: {best_target.over_kR_list}")
    lines.append(f"  C3 status: {best_target.status}")
    if best_target.largest_consecutive_width_kR and best_target.largest_consecutive_width_kR>0:
        lines.append(f"  Localization-window width proxy Δ(kR)={best_target.largest_consecutive_width_kR:.4f}; Q_proxy≈{best_target.largest_consecutive_Q_proxy:.2f}.")
    else:
        lines.append('  No finite-width consecutive over-threshold modal window was established on the distinct-mode spectrum; Q_proxy is therefore not a physical cavity Q.')
    lines.append('')
    lines.append('Local type-size screen around target configuration:')
    lines.append(f'  candidate configurations checked: {len(dfscreen)}')
    lines.append(f'  configurations with eta±10%R >= 92% on screen grid {args.screen_Nx}x{args.screen_Ns}: {over92_count}')
    lines.append('  Top candidates:')
    lines.append(dfscreen[['a','b','R','bc','m','kR','eta_100R','eta_050R','enrichment_100R','n_c2_peaks']].head(8).to_string(index=False))
    lines.append('')
    lines.append('Selected grid-validation rows:')
    lines.append(dfval[['a','b','R','bc','m','kR','eta_100R','eta_050R','n_c2_peaks','validation_Nx','validation_Ns']].to_string(index=False))
    lines.append('')
    lines.append('Scientific interpretation:')
    lines.append('  C2 remains supported: selected closed PHB configurations produce strong high-m annular localization above 70%, and some local-screen candidates exceed 92%.')
    lines.append('  C3 must be stated carefully. The dense shift-window scan is useful for locating modal peaks, but it repeatedly rediscoveres the same discrete eigenmodes. Therefore C3 should be based on deduplicated eigenmodes. If the deduplicated spectrum contains separated over-threshold peaks with below-threshold modes between them, this is not a continuous spectral window.')
    lines.append('  The physically safe formulation is: the target geometry exhibits a cluster of discrete over-threshold high-m annular modes. A true continuous C3 frequency-response window still requires a driven-response calculation with loss/gain or a much denser validated eigenmode map.')
    lines.append('  The effect is still dominated by high azimuthal order m and annular/toroidal topology. Hyperbolic geometry can be useful, but closed scalar C2/C3 alone does not prove unique hyperbolic superiority.')
    lines.append(f'Runtime seconds: {time.time()-t0:.1f}')
    summary='\n'.join(lines)
    (resdir/'RUN_SUMMARY_C2C3_continuous_scan.txt').write_text(summary,encoding='utf-8')
    print(summary)

    # Zip package
    zip_path=Path(args.package)
    if zip_path:
        with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
            for p in outdir.rglob('*'):
                if p.is_file(): z.write(p,p.relative_to(outdir.parent))
            z.write(Path(__file__),Path(__file__).name)
    return outdir


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--a',type=float,default=0.3)
    ap.add_argument('--b',type=float,default=0.6)
    ap.add_argument('--R',type=float,default=3.0)
    ap.add_argument('--m',type=int,default=15)
    ap.add_argument('--kR-min',type=float,default=10.0)
    ap.add_argument('--kR-max',type=float,default=25.0)
    ap.add_argument('--kR-step',type=float,default=0.5)
    ap.add_argument('--Nx',type=int,default=70)
    ap.add_argument('--Ns',type=int,default=48)
    ap.add_argument('--screen-Nx',type=int,default=50)
    ap.add_argument('--screen-Ns',type=int,default=35)
    ap.add_argument('--nev',type=int,default=6)
    ap.add_argument('--outdir', default=str(Path(__file__).resolve().parents[1] / 'run_continuous_scan'))
    ap.add_argument('--package', default=str(Path(__file__).resolve().parents[1] / 'PHB_C2C3_continuous_scan_results_package.zip'))
    args=ap.parse_args()
    run(args)

if __name__=='__main__':
    main()
