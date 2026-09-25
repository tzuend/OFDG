#!/usr/bin/env python3
"""Five-method errors, baseline preservation, and explicit metric definitions."""
import csv,json,math,collections
from pathlib import Path
from fe_mapping_reference import derivatives_at
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/reused-derivatives'
METHODS=['recursive','direct','fe_geometry','reused','reused_centered']
def read(p):return list(csv.DictReader(p.open()))
def write(name,rows):
 with (OUT/name).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def key(r):return tuple(r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative','method'])
def main():
 rows=read(OUT/'results.csv');lookup={key(r):r for r in rows};old=read(ROOT/'measurements/fe-mapping-derivatives/results.csv');delta=0
 for r in old:
  s=lookup[key(r)]
  for k in ['error_l2','relative_l2','sampled_linf','face_rms','jump_rms','input_error_l2','operator_error_l2']:
   a,b=float(r[k]),float(s[k]);delta=max(delta,abs(a-b)/max(1.,abs(a))) if math.isfinite(a) else delta
 assert delta==0,delta
 controls={k:max(float(r[k]) for r in rows) for k in ['reused_inverse_residual','reused_first_residual','reused_geometry_residual']}
 controls.update(baseline_rows=len(old),baseline_max_change=delta)
 controls['element_scaled_direct_agreement']=max(float(r['operator_error_l2'])/int(r['n'])**len(r['derivative']) for r in rows if r['method'].startswith('reused'))
 assert controls['element_scaled_direct_agreement']<1e-8
 assert controls['reused_inverse_residual']<1e-10 and controls['reused_first_residual']<1e-8 and controls['reused_geometry_residual']<1e-10
 basis=read(OUT/'basis-controls.csv');controls.update(basis_cases=len(basis),basis_polynomial_error=max(float(r['max_reference_error']) for r in basis),basis_constant_error=max(float(r['max_constant_error']) for r in basis),constant_coefficient_error=max(float(r['max_q_error']) for r in basis))
 for r in rows:
  order=len(r['derivative']);h=1/int(r['n']);r['element_scaled_l2']=h**order*float(r['error_l2']);r['element_scaled_operator_l2']=h**order*float(r['operator_error_l2'])
 write('errors.csv',rows)
 high=[];groups=collections.defaultdict(list)
 for r in read(OUT/'high-order.csv'):
  a=tuple(int(r[k]) for k in ['dx','dy','dz']);ref=0. if r['field']=='constant' else float(derivatives_at(int(r['dimension']),int(r['curved']),r['x'],r['y'],r['z'])[a]);v=float(r['value']);error=abs(v-ref)
  r.update(reference=ref,absolute_error=error,relative_error=error/abs(ref) if abs(ref)>1e-60 else math.nan,scaled_error=error/max(1.,abs(ref)),element_scaled_error=error/int(r['n'])**int(r['order']))
  high.append(r)
  if r['capacity']=='10':groups[tuple(r[k] for k in ['dimension','n','curved','field','method','order'])].append(r)
 controls['high_inverse_identity']=max(float(r['reused_inverse_residual']) for r in high);assert controls['high_inverse_identity']<1e-10
 write('high-order-errors.csv',high);summ=[]
 for key_,rs in groups.items():
  r=dict(zip(['dimension','n','curved','field','method','order'],key_));truth=sum(v['reference']**2 for v in rs);error=sum(v['absolute_error']**2 for v in rs)
  r.update(max_absolute_error=max(v['absolute_error'] for v in rs),max_scaled_error=max(v['scaled_error'] for v in rs),max_element_scaled_error=max(v['element_scaled_error'] for v in rs),discrete_relative_l2=math.sqrt(error/truth) if truth>1e-120 else math.nan);summ.append(r)
 write('high-order-summary.csv',summ)
 oldhigh=read(ROOT/'measurements/fe-mapping-derivatives/high-order/results.csv');hk=lambda r:tuple(r[k] for k in ['dimension','n','curved','capacity','point','field','dx','dy','dz','method']);hl={hk(r):r for r in high};hd=max(abs(float(r['value'])-float(hl[hk(r)]['value'])) for r in oldhigh);assert hd==0,hd;controls.update(high_baseline_rows=len(oldhigh),high_baseline_change=hd)
 probes=read(OUT/'stored-reference.csv');controls['stored_probe_count']=len(probes)//5;controls['stored_step_sensitivity']=max(float(r['step_sensitivity']) for r in probes)
 psummary=[]
 for dim in ['2','3']:
  for n in ['8','32']:
   for order in ['1','3','5']:
    for method in METHODS:
     rs=[r for r in probes if (r['dimension'],r['n'],r['order'],r['method'])==(dim,n,order,method)]
     if not rs:continue
     psummary.append(dict(dimension=dim,n=n,order=order,method=method,scaled_evaluation_error=max(float(r['scaled_error']) for r in rs),absolute_evaluation_error=max(float(r['absolute_error']) for r in rs),absolute_representation_error=max(float(r['representation_error']) for r in rs),absolute_total_error=max(float(r['total_error']) for r in rs)))
 write('stored-reference-summary.csv',psummary)
 # Lower-order outputs must be insensitive to the requested maximum order.
 capacity=0
 for r in high:
  if r['capacity']!='5':continue
  other=dict(r);other['capacity']='10';v=float(hl[hk(other)]['value']);capacity=max(capacity,abs(float(r['value'])-v)/max(1.,abs(v)))
 controls['capacity5_vs10']=capacity
 (OUT/'summary.json').write_text(json.dumps(dict(controls=controls,rows=len(rows),high_rows=len(high),high_summary=summ,stored_summary=psummary),indent=2))
 # Figures retain all five methods and mark zeros only by a plotting floor.
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 (OUT/'figures').mkdir(exist_ok=True);colors=['#b45c25','#497caa','#905aaa','#9c9f23','#087b72'];labels=['Recursive','Analytic direct','FE inverse series','Reused','Reused + centered']
 fig,axes=plt.subplots(2,2,figsize=(10,5.5),layout='constrained')
 for row,dim in enumerate(['2','3']):
  for col,n in enumerate(['8','32']):
   ax=axes[row,col]
   for method,color,label in zip(METHODS,colors,labels):
    rs=[r for r in summ if (r['dimension'],r['n'],r['curved'],r['field'],r['method'])==(dim,n,'1','quintic',method)];rs.sort(key=lambda r:int(r['order']));ax.semilogy([int(r['order']) for r in rs],[max(1e-18,r['max_scaled_error']) for r in rs],color=color,label=label,marker='.',linewidth=1.2)
   ax.set(title=f'{dim}D curved, h=1/{n}',xlabel='Derivative order',ylabel='max |error| / max(1, |truth|)');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=7,ncol=2);fig.savefig(OUT/'figures/high-orders.png',dpi=180);plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(10,3),layout='constrained')
 for ax,dim in zip(axes,['2','3']):
  for method,color,label in zip(METHODS,colors,labels):
   rs=[r for r in psummary if (r['dimension'],r['n'],r['method'])==(dim,'32',method)];ax.semilogy([int(r['order']) for r in rs],[max(1e-18,r['scaled_evaluation_error']) for r in rs],marker='o',color=color,label=label)
  ax.set(title=f'{dim}D, h=1/32: stored-function reference',xlabel='Derivative order',ylabel='max |evaluation error| / max(1, |reference|)');ax.grid(alpha=.2)
 axes[0].legend(fontsize=7);fig.savefig(OUT/'figures/stored-errors.png',dpi=180);plt.close(fig)
 print(json.dumps(controls,indent=2))
if __name__=='__main__':main()
