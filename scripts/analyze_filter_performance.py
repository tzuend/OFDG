#!/usr/bin/env python3
"""Summarize repeated production-kernel timings without changing study assets."""
import argparse,csv,json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/ofdg-performance-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
METHODS=['dg','ofdg','ofdg-kxrcf']
LABELS={'dg':'DG','ofdg':'OFDG','ofdg-kxrcf':'OFDG–KXRCF'}
COLORS={'dg':'#697586','ofdg':'#d97622','ofdg-kxrcf':'#007f86'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'measurements/performance/2026-09-07-main')
    parser.add_argument('--output',type=Path,default=ROOT/'output/performance')
    args=parser.parse_args(); source=args.input;out=args.output;out.mkdir(parents=True,exist_ok=True)
    state=json.loads((source/'results.json').read_text());summary=[];checks=[]
    med=lambda v:float(np.median(v))
    for run in state['runs']:
        if run['status']!='ok':continue
        samples=[s for s in run.get('samples',[]) if s['repeat']>=0]
        assert len(samples)==15,run['id']
        base={s['repeat']:s for s in samples if s['method']=='dg'}
        for method in METHODS:
            ss=[s for s in samples if s['method']==method];s=ss[0]
            times=np.array([s['total_seconds']/s['steps'] for s in ss])
            ratios=np.array([(s['total_seconds']/s['steps'])/(base[s['repeat']]['total_seconds']/base[s['repeat']]['steps']) for s in ss])
            row=dict(n=int(s['n']),p=int(s['order']),curved=bool(s['curved']),jump=bool(s['jump']),method=method,
                     cells=int(s['cells']),dofs=int(s['dofs']),time_ms=med(times)*1000,
                     time_q25_ms=float(np.quantile(times,.25))*1000,time_q75_ms=float(np.quantile(times,.75))*1000,
                     slowdown=med(ratios),ratio_q25=float(np.quantile(ratios,.25)),ratio_q75=float(np.quantile(ratios,.75)),
                     overhead_percent=(med(ratios)-1)*100,
                     direct_filter_share_percent=100*med([(s['detector_seconds']+s['decay_seconds'])/s['total_seconds'] for s in ss]),
                     dg_ms=1000*med([s['dg_seconds']/s['steps'] for s in ss]),
                     detector_ms=1000*med([s['detector_seconds']/s['steps'] for s in ss]),
                     decay_ms=1000*med([s['decay_seconds']/s['steps'] for s in ss]),
                     active_percent=100*med([s['active_count']/s['steps']/s['cells'] for s in ss]) if method=='ofdg-kxrcf' else None,
                     filter_setup_ms=1000*med([s['filter_setup'] for s in ss]),common_setup_ms=1000*s['common_setup'],
                     peak_mib=run['peak_mib'],repetitions=len(ss),trajectory_time=s['tf'],
                     steps_per_trajectory=s['steps']/s['cycles'])
            row['relative_iqr_percent']=100*(row['time_q75_ms']-row['time_q25_ms'])/row['time_ms']
            assert max(v['l2_error'] for v in ss)-min(v['l2_error'] for v in ss)<1e-12
            assert all(v['mass_drift']<1e-9 for v in ss)
            remainder=[v['total_seconds']-v['dg_seconds']-v['detector_seconds']-v['decay_seconds'] for v in ss]
            assert min(remainder)>=-1e-8
            row['other_ms']=1000*med([x/v['steps'] for x,v in zip(remainder,ss)])
            summary.append(row)
    def select(p=2,jump=False,curved=None,method=None,n=None):
        return sorted([r for r in summary if r['p']==p and r['jump']==jump and (curved is None or r['curved']==curved) and (method is None or r['method']==method) and (n is None or r['n']==n)],key=lambda r:r['n'])
    slopes=[]
    for p in [1,2,3]:
        for curved in [False,True]:
            for method in METHODS:
                rows=select(p=p,curved=curved,method=method)
                if len(rows)<2:continue
                exponent=float(np.polyfit(np.log([r['cells'] for r in rows]),np.log([r['time_ms'] for r in rows]),1)[0])
                local=[math.log(b['time_ms']/a['time_ms'])/math.log(b['cells']/a['cells']) for a,b in zip(rows,rows[1:])]
                slopes.append(dict(p=p,curved=curved,method=method,exponent=exponent,adjacent_exponents=local,n_range=[rows[0]['n'],rows[-1]['n']]))
    with (out/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,list(summary[0]));w.writeheader();w.writerows(summary)
    result=dict(summary=summary,slopes=slopes,failed=[r['id'] for r in state['runs'] if r['status']!='ok'],
                elapsed_seconds=state['elapsed_seconds'],prior_budget_seconds=state.get('prior_budget_seconds',0),
                sampled_peak_mib=max(r['peak_mib'] for r in state['runs']))
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18})
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(out/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9.3,3.7),sharey=True)
    for ax,curved in zip(axs,[False,True]):
        for method in METHODS[1:]:
            rows=select(curved=curved,method=method);v=np.array([r['slowdown'] for r in rows])
            ax.errorbar([r['cells'] for r in rows],v,yerr=[v-[r['ratio_q25'] for r in rows],[r['ratio_q75'] for r in rows]-v],
                        label=LABELS[method],color=COLORS[method],marker='o' if method=='ofdg' else 's',linestyle='-' if method=='ofdg' else '--',capsize=3)
        ax.axhline(1,color=COLORS['dg'],linestyle=':',label='DG baseline')
        ax.set_xscale('log',base=4); ns=sorted(set(r['cells'] for r in select(curved=curved)));ax.set_xticks(ns,labels=[str(n) for n in ns]);
        ax.set(title='Curved cubic geometry' if curved else 'Affine geometry',xlabel='Total cells',ylabel='Time per step / matched DG')
        ax.legend(fontsize=9)
    fig.suptitle('P2 smooth advection: median slowdown and interquartile range');fig.tight_layout();save(fig,'overhead')
    fig,axs=plt.subplots(1,2,figsize=(9.3,3.7))
    for ax,curved in zip(axs,[False,True]):
        for method in METHODS:
            rows=select(curved=curved,method=method)
            ax.loglog([r['cells'] for r in rows],[r['time_ms'] for r in rows],color=COLORS[method],marker='o',label=LABELS[method])
        rows=select(curved=curved,method='dg');xs=np.array([r['cells'] for r in rows]);ys=xs/xs[0]*rows[0]['time_ms']
        ax.plot(xs,ys,'k:',lw=1,label='Linear cell-count slope')
        ax.set(title='Curved' if curved else 'Affine',xlabel='Total cells',ylabel='Milliseconds / RK4 + filter step');ax.legend(fontsize=8.5)
    fig.suptitle('P2 smooth advection: empirical cost growth');fig.tight_layout();save(fig,'scaling')
    common=sorted(set(r['n'] for r in select(curved=False))&set(r['n'] for r in select(curved=True)))
    n=common[-1];fig,axs=plt.subplots(1,2,figsize=(9.3,3.8),sharey=True)
    for ax,curved in zip(axs,[False,True]):
        rows=[select(curved=curved,method=m,n=n)[0] for m in METHODS];bottom=np.zeros(3)
        for field,title,color in [('dg_ms','DG update','#697586'),('detector_ms','KXRCF','#007f86'),('decay_ms','OFDG decay','#d97622'),('other_ms','Loop / clock overhead','#ccd2d6')]:
            values=np.array([r[field] for r in rows]);ax.bar(range(3),values,bottom=bottom,label=title,color=color);bottom+=values
        ax.set_xticks(range(3),labels=[LABELS[m] for m in METHODS]);ax.set(title='Curved' if curved else 'Affine',ylabel='Milliseconds / step')
    handles,labels=axs[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.91),ncol=4,fontsize=8.5)
    maximum=max(sum(r[k] for k in ['dg_ms','detector_ms','decay_ms','other_ms']) for r in select(n=n))
    axs[0].set_ylim(0,maximum*1.10)
    fig.suptitle(f'P2 / {n} × {n} cells: directly measured time split')
    fig.tight_layout(rect=(0,0,1,.81));save(fig,'breakdown')
    lines=['# OFDG performance benchmark','',
           'Single-rank, release-build measurements on the current laptop. This is scalar 2D transport, not Euler or MPI scaling.', '',
           'Across the tested meshes, time per step grows approximately linearly with total cell count. For smooth P2 transport, OFDG adds about 44–50% over affine DG and 75–86% over curved DG. OFDG–KXRCF adds about 35–39% and 72–78%, respectively. These are stepping costs; curved operator construction is a separate, substantial initial cost.', '',
           'KXRCF is not uniformly faster: at 64² cells with curved P1 geometry it costs about 2.05 times DG, compared with 1.80 times DG for ungated OFDG. The following measurements separate degree, geometry, data dependence, and setup.', '',
           '## Method', '',
           '- Same production MFEM DG, RK4, OFDG decay, and KXRCF operators as the numerical drivers; production solver code is unchanged.',
           '- One warm-up round excluded from summaries, plus five measured repetitions. Method order rotates between repetitions. Each sample repeats a reset trajectory until at least 0.25 seconds of timed evolution has accumulated.',
           '- Matched physical interval T=0.005, CFL=0.15, velocity (0.7,0.3), post-step filtering. At fixed mesh and degree, every method uses the same timestep and number of steps. Short trajectories measure cost, not final solution accuracy.',
           '- Smooth initial data: sin²(pi(x+y)). Discontinuous initial data: 1 for x<0.5, otherwise 1.4, periodically extended. Curved geometry is the cubic interpolation of the existing 0.04*sin(2*pi*x)*sin(2*pi*y) displacement in both coordinates.',
           '- Constant prescribed velocity permits a fixed CFL timestep; the driver’s per-step CFL reduction is omitted equally for all methods. This is single-rank RK4-plus-filter cost rather than complete executable wall time.',
           '- Per-step clocks surround DG updates, indicator plus activity counting, and decay plus output copy. Setup, state reset, error checks, file I/O, process startup, and destruction are excluded from time-stepping measurements.',
           '- Slowdown is the median of five paired time-per-step ratios against DG. Error bars are interquartile ranges, not confidence intervals. Direct filter share is (indicator + decay) / total for that filtered run; it is different from overhead relative to a separate DG run.',
           '- Filter setup is the median repeated constructor time after warm-up; common mesh/DG setup is measured once per configuration. Setup is not included in stepping ratios. The process peak is sampled aggregate RSS.',
           '- Single MPI rank, one numerical-library thread. Eight available processors were confirmed; this uses less than half. Limits: 4 GB, 300 seconds/check, no new launches after 45 minutes, 60 minutes total. The main ledger reserves 300 seconds for the preliminary longer-trajectory pilot.', '',
           '## P2 mesh refinement, smooth data','',
           '| Geometry | Grid | DG ms/step | OFDG / DG | Gated / DG | OFDG share | Gated share | Active cells |',
           '|---|---:|---:|---:|---:|---:|---:|---:|']
    for curved in [False,True]:
        for n in sorted(set(r['n'] for r in select(curved=curved))):
            rs={m:select(curved=curved,method=m,n=n)[0] for m in METHODS}
            d,o,g=[rs[m] for m in METHODS]
            lines.append(f"| {'Curved' if curved else 'Affine'} | {n}² | {d['time_ms']:.3f} | {o['slowdown']:.3f}× | {g['slowdown']:.3f}× | {o['direct_filter_share_percent']:.1f}% | {g['direct_filter_share_percent']:.1f}% | {g['active_percent']:.2f}% |")
    lines+=['', '## Curved versus affine absolute stepping cost (P2)', '',
            'Each entry divides the curved median time per step by the affine median for the same method and grid. These ratios include the geometry cost in both the DG operator and filtering.', '',
            '| Grid | DG curved / affine | OFDG curved / affine | Gated curved / affine |',
            '|---|---:|---:|---:|']
    for n in sorted(set(r['n'] for r in select(curved=False)) & set(r['n'] for r in select(curved=True))):
        ratios=[select(curved=True,method=m,n=n)[0]['time_ms']/select(curved=False,method=m,n=n)[0]['time_ms'] for m in METHODS]
        lines.append(f"| {n}² | "+' | '.join(f'{v:.3f}×' for v in ratios)+' |')
    lines+=['','![Paired overhead](overhead.png)','','![Time split](breakdown.png)','','## Empirical growth','','Fit: time per step proportional to (total cells)^alpha. Alpha=1 is linear growth. These finite-range fits are not proofs of asymptotic complexity. At fixed final time and a CFL timestep, the number of steps also grows with cells per direction.','','| Degree | Geometry | DG alpha | OFDG alpha | Gated alpha |','|---|---|---:|---:|---:|']
    for p in [1,2,3]:
        for curved in [False,True]:
            ss=[next((s for s in slopes if s['p']==p and s['curved']==curved and s['method']==m),None) for m in METHODS]
            if all(ss):lines.append(f"| {p} | {'Curved' if curved else 'Affine'} | "+' | '.join(f"{s['exponent']:.3f}" for s in ss)+' |')
    lines+=['','![Scaling](scaling.png)','','## Degree dependence at 64² cells','','| Degree | Geometry | DG ms/step | OFDG / DG | Gated / DG |','|---|---|---:|---:|---:|']
    for p in [1,2,3]:
        for curved in [False,True]:
            rs=[select(p=p,curved=curved,method=m,n=64) for m in METHODS]
            if all(rs):lines.append(f"| {p} | {'Curved' if curved else 'Affine'} | {rs[0][0]['time_ms']:.3f} | {rs[1][0]['slowdown']:.3f}× | {rs[2][0]['slowdown']:.3f}× |")
    lines+=['','## KXRCF data dependence (P2)','','| Geometry | Grid | Data | Active cells | Gated / DG | Indicator ms | Decay ms |','|---|---:|---|---:|---:|---:|---:|']
    for curved in [False,True]:
        for n in [16,64]:
            for jump in [False,True]:
                rr=select(curved=curved,n=n,jump=jump,method='ofdg-kxrcf')
                if rr:
                    r=rr[0];lines.append(f"| {'Curved' if curved else 'Affine'} | {n}² | {'Jump' if jump else 'Smooth'} | {r['active_percent']:.2f}% | {r['slowdown']:.3f}× | {r['detector_ms']:.3f} | {r['decay_ms']:.3f} |")
    lines+=['','## Incremental setup (P2)','','| Geometry | Grid | Common setup, one observation | OFDG construction | OFDG + KXRCF construction |','|---|---:|---:|---:|---:|']
    for curved in [False,True]:
        for n in sorted(set(r['n'] for r in select(curved=curved))):
            d=select(curved=curved,method='dg',n=n)[0];o=select(curved=curved,method='ofdg',n=n)[0];g=select(curved=curved,method='ofdg-kxrcf',n=n)[0]
            lines.append(f"| {'Curved' if curved else 'Affine'} | {n}² | {d['common_setup_ms']:.1f} ms | {o['filter_setup_ms']:.1f} ms | {g['filter_setup_ms']:.1f} ms |")
    lines+=['','## Reproduce','',
            'From the repository root: `make build/release/benchmarks/filter_performance`, then `python3 scripts/run_filter_performance.py --output measurements/performance/new-run --prior-budget-seconds 0`. Finally run `python3 scripts/analyze_filter_performance.py --input measurements/performance/new-run`. Use a new directory when sources, settings, or binaries change. The current run used a 300-second conservative pilot reservation.', '',
            '## Interpretation and limits','',
            'KXRCF can skip inactive face contributions and shell updates, but the current OFDG kernel still computes global scaling and builds derivative states on every cell. Zero activated cells therefore does not mean zero filtering cost. Curved geometry additionally increases physical quadrature/evaluation work. These code observations explain plausible costs; the present phase timers do not resolve every internal sub-operation.', '',
            'Overhead relative to DG and direct filtering share are different measurements. The DG phase inside a filtered run may itself take longer than in the separate baseline, so the slowdown need not equal one divided by the unfiltered time fraction. Cache effects or processor variation are possible contributors, but these measurements do not distinguish them. The largest time-per-step interquartile range is about 14% of its median; small differences should not be overinterpreted.', '',
            'Comparisons are specific to scalar quadrilateral DG, these degrees, cubic deformation, and one rank. They do not establish Euler overhead, performance on triangles/3D, strong scaling, or large-system asymptotic behaviour. The ordinary DG update itself can have a different curved-mesh cost; compare both absolute times and matched DG ratios. Repeated short reset trajectories emphasize early-time detector activity.', '',
            f"Main benchmark ledger: {state['elapsed_seconds']/60:.2f} minutes including {state.get('prior_budget_seconds',0)/60:.1f} minutes reserved for pilot work. Sampled peak RSS: {result['sampled_peak_mib']:.1f} MiB. Completed configurations: {sum(r['status']=='ok' for r in state['runs'])}/{len(state['configs'])}. Failed/limited configurations: {result['failed'] or 'none'}.", '',
            f"Raw data and exact commands: `{source.resolve()}`. Summary: `summary.csv`; fitted exponents and quartiles: `summary.json`. The pilot remains in the sibling `2026-09-07` directory. No thesis result asset was replaced."]
    (out/'benchmark.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(configurations=len(state['runs']),rows=len(summary),slopes=slopes,elapsed=state['elapsed_seconds']),indent=2))
if __name__=='__main__':main()
