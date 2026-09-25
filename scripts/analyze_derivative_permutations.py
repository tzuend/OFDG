#!/usr/bin/env python3
"""Analyze the isolated study; write derived tables and readable difference plots."""
import csv,json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/ofdg-derivative-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/derivative-permutations'
def main():
    rows=list(csv.DictReader((OUT/'results.csv').open()))
    for r in rows:
        for k in r.keys()-{'family','field','multi','path'}:r[k]=float(r[k])
    def key(r):return tuple(r[k] for k in ['family','n','p','amplitude','qextra','field','multi'])
    groups={}
    for r in rows:groups.setdefault(key(r),{})[r['path']]=r
    comparisons=[]
    for keyval,g in groups.items():
        c=g[keyval[-1]];a=g['mean'];e=c['error_l2']
        comparisons.append({**{k:c[k] for k in ['family','n','p','amplitude','qextra','field','multi']},'canonical_error':e,'mean_error':a['error_l2'],'mean_to_canonical':a['error_l2']/e if e else 1.,'change_over_error':c['path_minus_mean_l2']/e if e else 0.,'change_l2':c['path_minus_mean_l2'],'canonical_relative':e/c['exact_l2'],'floor':c['best_projection_error_l2']})
    with (OUT/'comparisons.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=comparisons[0]);w.writeheader();w.writerows(comparisons)
    fine=[r for r in comparisons if r['n']==8 and r['amplitude']>0 and r['qextra']==0 and len(r['multi'])==3]
    affine=[r for r in comparisons if r['amplitude']==0]
    summary={'runs':json.loads((OUT/'run.json').read_text()),'comparison_count':len(comparisons),'fine_third_count':len(fine),
       'fine_better':sum(r['mean_to_canonical']<1 for r in fine),'fine_ratio_min':min(r['mean_to_canonical'] for r in fine),'fine_ratio_median':float(np.median([r['mean_to_canonical'] for r in fine])),'fine_ratio_max':max(r['mean_to_canonical'] for r in fine),
       'fine_change_min':min(r['change_over_error'] for r in fine),'fine_change_max':max(r['change_over_error'] for r in fine),
       'affine_max_change':max(r['change_l2'] for r in affine),'affine_cubic_error':max(r['canonical_error'] for r in affine if r['field']=='cubic'),
       'first_weak_max':max(r['first_weak_relative'] for r in rows),'constant_max':max(r['constant_max'] for r in rows),
       'geometry_min_global_det':min(r['min_jacobian_det']*r['n']**2 for r in rows)}
    lookup={key(r)+(r['path'],):r for r in rows};qchanges=[]
    for r in rows:
        if r['qextra']==8:
            k=list(key(r));k[4]=0;b=lookup[tuple(k)+(r['path'],)]
            qchanges.append(abs(r['error_l2']-b['error_l2'])/max(b['error_l2'],1e-12))
    summary['quadrature_max_relative_error_change']=max(qchanges)
    summary['representative']=[r for r in comparisons if r['n']==8 and r['amplitude']==1.2 and r['qextra']==0 and r['multi']=='xxy' and r['field']=='sine']
    summary['refinement']=[]
    for fam in ['quad','tri']:
        for p in [3,4,5]:
            series=sorted([r for r in comparisons if r['family']==fam and r['p']==p and r['amplitude']==1.2 and r['qextra']==0 and r['multi']=='xxy' and r['field']=='sine'],key=lambda r:r['n'])
            a,b=series[-2:];summary['refinement'].append({'family':fam,'p':p,'canonical_rate':math.log2(a['canonical_error']/b['canonical_error']),'mean_rate':math.log2(a['mean_error']/b['mean_error'])})
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    figures=OUT/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(10,5.7),layout='constrained')
    for i,fam in enumerate(['quad','tri']):
        for j,p in enumerate([3,4,5]):
            ax=axes[i,j];data=sorted([r for r in comparisons if r['family']==fam and r['p']==p and r['amplitude']==1.2 and r['qextra']==0 and r['multi']=='xxy' and r['field']=='sine'],key=lambda r:r['n'])
            ns=[r['n'] for r in data]
            ax.loglog(ns,[r['canonical_error'] for r in data],'o-',color='#007f86',label='xxy')
            ax.loglog(ns,[r['mean_error'] for r in data],'x--',color='#c05a14',label='permutation mean')
            ax.set_title(f'{"Quadrilateral Q" if fam=="quad" else "Triangle P"}{p}');ax.set_xticks(ns,labels=[str(int(n)) for n in ns]);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());ax.grid(alpha=.2,which='both');ax.set_xlabel('n (mesh spacing 1/n)')
            if j==0:ax.set_ylabel('Physical L2 error')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=2)
    fig.savefig(figures/'convergence.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.3),layout='constrained')
    for ax,fam in zip(axes,['quad','tri']):
        ax.axhline(1,color='#555',lw=1)
        for field,marker,col in [('sine','o','#007f86'),('exponential','^','#c05a14'),('cubic','s','#7351a0')]:
            rr=[r for r in fine if r['family']==fam and r['field']==field]
            ax.scatter([100*r['change_over_error'] for r in rr],[r['mean_to_canonical'] for r in rr],marker=marker,c=col,s=32,alpha=.75,label=field)
        ax.set_xlabel('Change from xxy / xyy to mean\n(% of canonical derivative error)');ax.set_title('Quadrilaterals' if fam=='quad' else 'Triangles');ax.grid(alpha=.15)
    axes[0].set_ylabel('Mean error / canonical error\n(< 1: averaging helps)')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3)
    fig.savefig(figures/'averaging.png',dpi=180);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['runs']},indent=2))
if __name__=='__main__':main()
