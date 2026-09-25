#!/usr/bin/env python3
"""Analyse analytic derivative errors; no numerical solutions are modified."""
import csv,json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/ofdg-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/derivative-accuracy'
def main():
 rows=list(csv.DictReader((OUT/'results.csv').open()))
 textkeys={'family','basis','map','field','derivative'}
 for r in rows:
  for k in r:
   if k not in textkeys:r[k]=float(r[k])
  r['order']=len(r['derivative'])
 def pick(**kwargs):return [r for r in rows if all(r[k]==v for k,v in kwargs.items())]
 base=[r for r in rows if r['qextra']==0 and r['evalextra']==0]
 summary={'rows':len(rows),'run':json.loads((OUT/'run.json').read_text()),'weak_max':max(r['weak_residual'] for r in rows),'constant_linf_max':max(r['sampled_linf'] for r in rows if r['field']=='constant'),'represented_input_max':max([r['input_error_l2'] for r in rows if r['field']=='represented_quadratic'] or [0]),'representatives':[],'isolated':[]}
 # Error envelopes use maxima over all canonical derivatives at the indicated order.
 for mapping,fields,amplitude,key in [('coupled',['sine','shortwave','exponential'],1.2,'representatives'),('separable',['represented_quadratic'],.7,'isolated')]:
  for family in ['quad','tri']:
   for p in [2,3,4]:
    for field in fields:
     for order in range(1,min(p,3)+1):
      rr=[r for r in base if r['map']==mapping and r['family']==family and r['p']==p and r['field']==field and r['order']==order and r['n']==8 and r['basis']=='gll' and r['amplitude']==amplitude]
      if not rr:continue
      worst=max(rr,key=lambda r:r['relative_l2']);summary[key].append({'family':family,'p':p,'field':field,'order':order,'derivative':worst['derivative'],'abs_l2':worst['error_l2'],'rel_l2':worst['relative_l2'],'digits':-math.log10(worst['relative_l2']),'face_rel_max':max(r['face_relative_l2'] for r in rr),'sample_rel_max':max(r['relative_sampled_linf'] for r in rr)})
 # Controls compare identical geometry/order/field/derivative; report absolute changes too.
 lookup={(r['family'],r['n'],r['p'],r['amplitude'],r['map'],r['field'],r['derivative']):r for r in base if r['basis']=='gll'}
 comparisons=[]
 for r in rows:
  key=(r['family'],r['n'],r['p'],r['amplitude'],r['map'],r['field'],r['derivative']);ref=lookup.get(key)
  if ref is None:continue
  kind='basis' if r['basis']!='gll' else ('assembly' if r['qextra'] else ('evaluation' if r['evalextra'] else None))
  if not kind:continue
  change=abs(r['error_l2']-ref['error_l2']);comparisons.append({'kind':kind,'family':r['family'],'map':r['map'],'p':r['p'],'field':r['field'],'derivative':r['derivative'],'basis':r['basis'],'absolute_error_norm_change':change,'scaled_change':change/max(ref['truth_l2'],ref['input_l2'],1e-30),'relative_error_norm_change':change/ref['error_l2'] if ref['error_l2']>1e-10 else float('nan')})
 summary['controls']={k:{'max_scaled_change':max([c['scaled_change'] for c in comparisons if c['kind']==k] or [0]),'max_absolute_change':max([c['absolute_error_norm_change'] for c in comparisons if c['kind']==k] or [0]),'max_nonconstant_scaled_change':max([c['scaled_change'] for c in comparisons if c['kind']==k and c['field']!='constant'] or [0]),'max_relative_error_change':max([c['relative_error_norm_change'] for c in comparisons if c['kind']==k and math.isfinite(c['relative_error_norm_change'])] or [0])} for k in ['basis','assembly','evaluation']}
 # Orthogonality decomposition and exact-input first-derivative excess.
 summary['pythagorean_relative_defect']=max(abs(r['error_l2']**2-r['best_projection_l2']**2-r['excess_over_projection_l2']**2)/r['error_l2']**2 for r in rows if r['field']!='constant' and r['error_l2']>1e-9)
 summary['isolated_first_excess_max']=max([r['excess_over_projection_l2'] for r in rows if r['field']=='represented_quadratic' and r['order']==1] or [0])
 rates=[]
 for r in base:
  if r['basis']!='gll' or r['field']=='constant' or not r['amplitude']:continue
  prev=lookup.get((r['family'],r['n']/2,r['p'],r['amplitude'],r['map'],r['field'],r['derivative']))
  if prev and min(prev['error_l2'],r['error_l2'])>1e-10:rates.append({'family':r['family'],'map':r['map'],'field':r['field'],'p':r['p'],'n':r['n'],'amplitude':r['amplitude'],'derivative':r['derivative'],'rate':math.log2(prev['error_l2']/r['error_l2'])})
 summary['rates']=rates
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
 for name,data in [('controls',comparisons),('convergence',rates),('representative_errors',summary['representatives']+summary['isolated'])]:
  if data:
   with (OUT/f'{name}.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
 colours={1:'#777777',2:'#007f86',3:'#b65b24',4:'#6c4da0'}
 for mapping,field,amplitude,name in [('coupled','sine',1.2,'accuracy'),('separable','represented_quadratic',.7,'isolated')]:
  fig,axes=plt.subplots(2,3,figsize=(10,5.4),sharex='col')
  for row,family in enumerate(['quad','tri']):
   for col,order in enumerate([1,2,3]):
    ax=axes[row,col]
    for p in [2,3,4]:
     if p<order:continue
     selected=[r for r in base if r['family']==family and r['map']==mapping and r['field']==field and r['amplitude']==amplitude and r['basis']=='gll' and r['p']==p and r['order']==order]
     ns=sorted(set(r['n'] for r in selected));vals=[max(r['relative_l2'] for r in selected if r['n']==n) for n in ns]
     if ns:ax.loglog(ns,vals,'o-',color=colours[p],label=('Q' if family=='quad' else 'P')+str(p),markersize=4)
    ax.set_title(('Quadrilaterals' if family=='quad' else 'Triangles')+f': order {order}');ax.grid(alpha=.2,which='both');ax.set_xticks([2,4,8,16],[2,4,8,16]);ax.legend(fontsize=8)
    if col==0:ax.set_ylabel('Largest relative L2 error')
    if row==1:ax.set_xlabel('Subdivisions per side, n')
  fig.tight_layout();fig.savefig(figdir/f'{name}.png',bbox_inches='tight');plt.close(fig)
 print(json.dumps({k:v for k,v in summary.items() if k not in ['run','rates','representatives','isolated']},indent=2))
if __name__=='__main__':main()
