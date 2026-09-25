#!/usr/bin/env python3
"""Matched-input physical reconstruction analysis; retains previous studies."""
import csv,json,math,os,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/ofdg-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/physical-reconstruction'
METHODS=['recursive','physical_k','physical_k1','physical_k2','oracle_k2']
def read(path):
 rows=list(csv.DictReader(path.open()))
 for r in rows:
  for k in r:
   if k not in ['family','basis','map','field','derivative','method']:r[k]=float(r[k])
  r['order']=len(r['derivative'])
 return rows

def main():
 rows=read(OUT/'results.csv');base=[r for r in rows if r['extra']==0 and r['basis']=='gll']
 def key(r):return tuple(r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative'])
 lookup={(key(r),r['method']):r for r in rows}
 comparisons=[]
 for r in rows:
  if r['method'] not in ['physical_k','physical_k1','physical_k2']:continue
  old=lookup[(key(r),'recursive')]
  if r['field']=='constant' or old['error_l2']<1e-9 or old['truth_l2']<1e-10:continue
  comparisons.append({k:r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative','order','method']}|{'error_ratio':r['error_l2']/old['error_l2'],'face_ratio':r['face_rms']/old['face_rms'] if old['face_rms']>1e-9 else float('nan'),'jump_ratio':r['jump_rms']/old['jump_rms'] if old['jump_rms']>1e-9 else float('nan'),'relative_l2':r['relative_l2'],'error_l2':r['error_l2']})
 reps=[]
 for mapping,field,amp in [('coupled','sine',1.2),('separable','represented_quadratic',.7)]:
  for family in ['quad','tri']:
   for p in [2,3,4]:
    for order in range(1,min(p,3)+1):
     for method in METHODS:
      rr=[r for r in base if r['map']==mapping and r['field']==field and r['amplitude']==amp and r['family']==family and r['p']==p and r['n']==8 and r['order']==order and r['method']==method]
      if rr:
       w=max(rr,key=lambda r:r['relative_l2']);reps.append({k:w[k] for k in ['map','field','family','p','order','method','derivative','relative_l2','error_l2']}|{'face_relative_max':max(r['face_relative_l2'] for r in rr)})
 aggregate=[]
 for mapping in ['coupled','separable']:
  for method in METHODS[1:4]:
   rr=[r for r in comparisons if r['map']==mapping and r['method']==method and r['n']==8 and r['basis']=='gll' and r['extra']==0 and r['amplitude']==(1.2 if mapping=='coupled' else .7)]
   if rr:
    aggregate.append({'map':mapping,'method':method,'count':len(rr),**{metric+suffix:float(func([r[metric] for r in rr if math.isfinite(r[metric])])) for metric in ['error_ratio','face_ratio','jump_ratio'] for suffix,func in [('_min',min),('_median',np.median),('_max',max)]},'improved':sum(r['error_ratio']<1-1e-6 for r in rr),'worsened':sum(r['error_ratio']>1+1e-6 for r in rr)})
 controls={}
 controls['mean_defect_max']=max(r['mean_defect'] for r in rows)
 controls['solve_residual_max']=max(r['solve_residual'] for r in rows)
 controls['constant_max_sampled']=max(r['sampled_linf'] for r in rows if r['field']=='constant')
 controls['affine_quadratic_max_l2']=max([r['error_l2'] for r in rows if r['field']=='physical_quadratic'] or [0])
 changes={'basis':[],'quadrature':[],'affine_identity':[],'old_baseline':[],'oracle_exact_input':[]}
 for r in rows:
  if r['basis']!='gll' or r['extra']:
   target=r.copy();target['basis']='gll';target['extra']=0.;ref=lookup.get((key(target),r['method']))
   if ref:changes['basis' if r['basis']!='gll' else 'quadrature'].append(abs(r['error_l2']-ref['error_l2'])/max(1.,ref['truth_l2']))
  if r['family']=='tri' and r['amplitude']==0 and r['method'].startswith('physical_'):
   ref=lookup[(key(r),'recursive')];changes['affine_identity'].append(abs(r['error_l2']-ref['error_l2'])/max(1.,ref['truth_l2']))
  if r['field']=='represented_quadratic' and r['method']=='physical_k2':
   ref=lookup[(key(r),'oracle_k2')];changes['oracle_exact_input'].append(abs(r['error_l2']-ref['error_l2'])/max(1.,ref['truth_l2']))
 previous=read(ROOT/'measurements/derivative-accuracy/results.csv')
 previous_lookup={tuple(r[k] for k in ['family','n','p','amplitude','basis','map','field','derivative']):r for r in previous if r['qextra']==0 and r['evalextra']==0}
 for r in rows:
  if r['method']!='recursive' or r['extra']:continue
  ref=previous_lookup.get(tuple(r[k] for k in ['family','n','p','amplitude','basis','map','field','derivative']))
  if ref:changes['old_baseline'].append(abs(r['error_l2']-ref['error_l2'])/max(1.,ref['truth_l2']))
 for name,values in changes.items():controls[name+'_max_scaled_norm_change']=max(values or [0]);controls[name+'_count']=len(values)
 # Nested physical spaces must not increase the solution reconstruction error.
 nesting=[]
 for r in rows:
  if r['method']!='physical_k':continue
  values=[lookup[(key(r),m)]['reconstruction_l2'] for m in METHODS[1:4]]
  nesting.extend([max(0.,values[i+1]-values[i]) for i in [0,1]])
 controls['nested_reconstruction_violation']=max(nesting or [0])
 assert controls['mean_defect_max']<1e-9
 assert controls['solve_residual_max']<1e-10
 assert controls['affine_quadratic_max_l2']<1e-7
 assert controls['affine_identity_max_scaled_norm_change']<1e-7
 assert controls['old_baseline_max_scaled_norm_change']<1e-7
 assert controls['nested_reconstruction_violation']<1e-10
 summary={'rows':len(rows),'run':json.loads((OUT/'run.json').read_text()),'representatives':reps,'aggregate':aggregate,'controls':controls}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
 for name,data in [('comparisons',comparisons),('representatives',reps),('aggregate',aggregate)]:
  if data:
   with (OUT/f'{name}.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':180})
 # Ratios rather than overlapping solution/error curves.
 fig,axes=plt.subplots(2,2,figsize=(10,5.6))
 for row,mapping in enumerate(['coupled','separable']):
  for col,family in enumerate(['quad','tri']):
   ax=axes[row,col];matrix=np.full((3,3),np.nan)
   for order in [1,2,3]:
    ref=next((r for r in reps if r['map']==mapping and r['family']==family and r['p']==3 and r['order']==order and r['method']=='recursive'),None)
    for j,method in enumerate(METHODS[1:4]):
     cur=next((r for r in reps if r['map']==mapping and r['family']==family and r['p']==3 and r['order']==order and r['method']==method),None)
     if ref and cur:matrix[order-1,j]=cur['relative_l2']/ref['relative_l2']
   im=ax.imshow(np.log10(matrix),cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
   for i in range(3):
    for j in range(3):
     if np.isfinite(matrix[i,j]):ax.text(j,i,f'{matrix[i,j]:.3g}x',ha='center',va='center',color='white' if abs(np.log10(matrix[i,j]))>1.2 else 'black')
   ax.set_xticks([0,1,2],['q=k','q=k+1','q=k+2']);ax.set_yticks([0,1,2],['1st','2nd','3rd']);ax.set_ylabel('Derivative order');ax.set_title(('Q3' if family=='quad' else 'P3')+' / '+('projected sine' if mapping=='coupled' else 'exactly represented input'))
 fig.suptitle('Error ratio to current recursion, n=8 (below 1 is better)',fontsize=12)
 fig.tight_layout(rect=(0,0,.91,.94));cax=fig.add_axes([.93,.15,.015,.65]);cb=fig.colorbar(im,cax=cax);cb.set_label('log10(error ratio)');fig.savefig(figdir/'ratios.png',bbox_inches='tight');plt.close(fig)
 # Directly plot improvements/deterioration in highest derivative under refinement.
 fig,axes=plt.subplots(1,2,figsize=(10,3.1))
 for ax,mapping in zip(axes,['coupled','separable']):
  for fam,marker in [('quad','o'),('tri','s')]:
   for method,colour in [('physical_k','#a35a2a'),('physical_k1','#007f86'),('physical_k2','#7653a0')]:
    rr=[r for r in comparisons if r['family']==fam and r['map']==mapping and r['p']==3 and r['derivative']=='xxx' and r['field']==('sine' if mapping=='coupled' else 'represented_quadratic') and r['basis']=='gll' and r['extra']==0 and r['amplitude']==(1.2 if mapping=='coupled' else .7) and r['method']==method]
    rr.sort(key=lambda r:r['n'])
    if rr:ax.plot([r['n'] for r in rr],[r['error_ratio'] for r in rr],marker=marker,color=colour,ls='-' if fam=='quad' else '--',label=('Q3' if fam=='quad' else 'P3')+' / '+method.replace('physical_','q=').replace('k1','k+1').replace('k2','k+2'))
  ax.axhline(1,color='black',lw=.8);ax.set_yscale('log');ax.set_xscale('log',base=2);ax.set_xticks([2,4,8,16],[2,4,8,16]);ax.set_xlabel('Subdivisions per side');ax.set_title('Projected sine' if mapping=='coupled' else 'Exactly represented input');ax.grid(alpha=.2);ax.set_ylabel('xxx error / current xxx error')
 axes[1].legend(fontsize=7,ncol=2);fig.tight_layout();fig.savefig(figdir/'refinement.png',bbox_inches='tight');plt.close(fig)
 print(json.dumps({'rows':len(rows),'controls':controls,'aggregate':aggregate},indent=2))
if __name__=='__main__':main()
