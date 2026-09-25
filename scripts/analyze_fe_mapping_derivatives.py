#!/usr/bin/env python3
"""Matched baseline comparisons, independent controls, and high-order limits."""
import csv,json,math,statistics,time
from pathlib import Path
from decimal import localcontext
from analyze_direct_derivatives import reference_derivative
from fe_mapping_reference import derivatives_at,finite_difference
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/fe-mapping-derivatives'
CATS={'family','basis','map','field','derivative','method'}
def read(path,cats=CATS):
 rows=list(csv.DictReader(path.open()))
 for r in rows:
  for k in r:
   if k not in cats:r[k]=float(r[k])
 return rows
def write(name,rows):
 with (OUT/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def key(r):return tuple(r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative','method'])
def main():
 start=time.monotonic();rows=read(OUT/'results.csv');lookup={key(r):r for r in rows};assert len(rows)==len(lookup)
 previous=read(ROOT/'measurements/direct-derivatives/results.csv');baseline=[]
 for r in previous:
  b=lookup[key(r)]
  for metric in ['error_l2','face_rms','jump_rms','operator_error_l2','input_error_l2']:
   baseline.append(abs(r[metric]-b[metric])/max(1.,r['truth_l2']))
 controls={k:max(r[k] for r in rows) for k in ['geometry_residual','geometry_derivative_residual','geometry_jacobian_residual','face_geometry_residual','fe_inverse_residual','method_agreement','fe_first_residual','representation_residual']}
 controls.update(baseline_rows=len(previous),baseline_max_scaled_difference=max(baseline))
 assert controls['baseline_max_scaled_difference']<1e-10
 for k in ['geometry_residual','geometry_derivative_residual','geometry_jacobian_residual','face_geometry_residual','fe_inverse_residual']:assert controls[k]<1e-10,(k,controls[k])
 for k in ['method_agreement','fe_first_residual']:assert controls[k]<1e-8,(k,controls[k])
 comparisons=[]
 for r in rows:
  if r['method']!='fe_geometry':continue
  q=dict(r);q['method']='direct';d=lookup[key(q)];q['method']='recursive';b=lookup[key(q)]
  item={k:r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative']}
  for metric in ['error_l2','face_rms','jump_rms']:
   item[metric+'_fe_over_recursive']=r[metric]/b[metric] if b[metric]>1e-9 else math.nan
   item[metric+'_fe_minus_direct']=r[metric]-d[metric]
  item['fe_operator_relative_l2']=r['operator_relative_l2'];comparisons.append(item)
 write('comparisons',comparisons)
 for name,selector in [('basis',lambda r:r['basis']!='gll'),('quadrature',lambda r:r['extra']>0)]:
  differences=[]
  for r in rows:
   if not selector(r):continue
   q=dict(r);q.update(basis='gll',extra=0.)
   b=lookup[key(q)];differences.append(abs(r['error_l2']-b['error_l2'])/max(1.,r['truth_l2']))
  controls[name+'_max_scaled_norm_change']=max(differences)
 controls['constant_fe_max_l2']=max(r['error_l2'] for r in rows if r['field']=='constant' and r['method']=='fe_geometry')
 controls['affine_quadratic_fe_max_l2']=max(r['error_l2'] for r in rows if r['field']=='physical_quadratic' and r['method']=='fe_geometry')
 assert controls['constant_fe_max_l2']<1e-6 and controls['affine_quadratic_fe_max_l2']<1e-6
 independent=[]
 for m,x,y,v in csv.reader((OUT/'control-probes.csv').open()):
  coarse=reference_derivative(m,int(x),int(y),'1e-5');fine=reference_derivative(m,int(x),int(y),'1e-6')
  independent.append(dict(dimension=2,map=m,derivative='x'*int(x)+'y'*int(y),value=float(v),reference=fine,scaled_error=abs(float(v)-fine)/max(1,abs(fine)),step_sensitivity=abs(coarse-fine)/max(1,abs(fine))))
 high_raw=list(csv.DictReader((OUT/'high-order/results.csv').open()));high=[];precision=[]
 for r in high_raw:
  dim=int(r['dimension']);curved=int(r['curved']);alpha=tuple(int(r[k]) for k in ['dx','dy','dz']);value=float(r['value'])
  ref=0. if r['field']=='constant' else float(derivatives_at(dim,curved,r['x'],r['y'],r['z'])[alpha])
  if abs(ref)<1e-60:ref=0.
  error=abs(value-ref);rr={k:(r[k] if k in {'method','field'} else float(r[k])) for k in r};rr.update(reference=ref,absolute_error=error,scaled_error=error/max(1,abs(ref)),relative_error=error/abs(ref) if ref else math.nan)
  high.append(rr)
  if r['field']=='quintic' and r['method']=='fe_geometry' and r['n']=='2' and r['curved']=='1' and r['capacity']=='10' and r['point']=='0':
   if int(r['order'])<=3:
    coarse=float(finite_difference(dim,curved,(r['x'],r['y'],r['z']),alpha,'1e-5'));fine=float(finite_difference(dim,curved,(r['x'],r['y'],r['z']),alpha,'1e-6'))
    independent.append(dict(dimension=dim,map='separable' if dim==2 else 'shear',derivative=''.join('xyz'[d]*alpha[d] for d in range(dim)),value=value,reference=fine,scaled_error=abs(value-fine)/max(1,abs(fine)),step_sensitivity=abs(coarse-fine)/max(1,abs(fine))))
   if int(r['order']) in [5,10]:
    finer=float(derivatives_at(dim,curved,r['x'],r['y'],r['z'],120)[alpha]);precision.append(abs(ref-finer)/max(1,abs(ref)))
 controls['independent_probe_count']=len(independent);controls['independent_max_scaled_error']=max(r['scaled_error'] for r in independent);controls['independent_max_step_sensitivity']=max(r['step_sensitivity'] for r in independent);controls['high_reference_precision_change']=max(precision)
 assert controls['independent_max_scaled_error']<1e-7 and controls['independent_max_step_sensitivity']<1e-7
 write('independent-controls',independent);write('high-order-errors',high)
 groups={}
 for r in high:
  if r['capacity']!=10:continue
  k=tuple(r[x] for x in ['dimension','n','curved','field','order','method']);groups.setdefault(k,[]).append(r)
 aggregate=[]
 for k,rs in groups.items():
  a=dict(zip(['dimension','n','curved','field','order','method'],k));a.update(max_scaled_error=max(r['scaled_error'] for r in rs),max_absolute_error=max(r['absolute_error'] for r in rs),max_nonzero_relative_error=max([r['relative_error'] for r in rs if math.isfinite(r['relative_error'])] or [math.nan]));aggregate.append(a)
 write('high-order-summary',aggregate)
 hk=lambda r:tuple(r[x] for x in ['dimension','n','curved','point','field','dx','dy','dz','method'])
 high10={hk(r):r for r in high if r['capacity']==10};overlap=[]
 for r in high:
  if r['capacity']==5:
   b=high10[hk(r)];overlap.append(abs(r['value']-b['value'])/max(1,abs(b['value'])))
 controls['order5_vs_order10_max_scaled_difference']=max(overlap)
 assert controls['order5_vs_order10_max_scaled_difference']<1e-8
 bymethod={}
 for dim in ['quad','tri','hex']:
  for field in (['represented_polynomial','sine'] if dim=='hex' else ['represented_quadratic','sine']):
   for method in ['recursive','direct','fe_geometry']:
    candidates=[r for r in rows if r['family']==dim and r['field']==field and r['method']==method and r['n']==(4 if dim=='hex' else 8) and r['p']==3 and r['basis']=='gll' and r['extra']==0 and r['amplitude']==(.2 if dim=='hex' else .7 if field=='represented_quadratic' else 1.2) and len(r['derivative'])==3]
    bymethod[f'{dim}-{field}-{method}']=max(r['relative_l2'] for r in candidates if math.isfinite(r['relative_l2']))
 summary=dict(controls=controls,row_count=len(rows),high_order_row_count=len(high),representative_worst_third_relative=bymethod,high_order_summary=aggregate,analysis_seconds=time.monotonic()-start)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'controls':controls,'representative_worst_third_relative':bymethod,'analysis_seconds':summary['analysis_seconds']},indent=2))
 plots(rows,aggregate)
def plots(rows,high):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import numpy as np
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
 methods=[('recursive','#9c4b28','o'),('direct','#5a67a1','s'),('fe_geometry','#00878a','^')]
 fig,axes=plt.subplots(1,2,figsize=(9,3.2),layout='constrained')
 for ax,field in zip(axes,['represented_quadratic','sine']):
  for method,color,marker in methods:
   ys=[]
   for n in [4,8,16]:
    rs=[r for r in rows if r['family']=='quad' and r['field']==field and r['n']==n and r['p']==3 and r['method']==method and r['basis']=='gll' and r['extra']==0 and r['amplitude']==(.7 if field=='represented_quadratic' else 1.2) and len(r['derivative'])==3];ys.append(max(r['relative_l2'] for r in rs))
   ax.plot([4,8,16],ys,color=color,marker=marker,label={'recursive':'Recursive','direct':'Analytic direct','fe_geometry':'FE geometry'}[method],linestyle='--' if method=='direct' else '-')
  ax.set_yscale('log');ax.set_xticks([4,8,16]);ax.set_xlabel('Subdivisions per direction');ax.set_title('Exactly represented quadratic' if field=='represented_quadratic' else 'Projected sine');ax.grid(alpha=.2);ax.legend(fontsize=8)
 axes[0].set_ylabel('Worst third-derivative relative L2 error');fig.savefig(figdir/'comparison-2d.png',dpi=180);plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(9,3.1),layout='constrained')
 for ax,dim in zip(axes,[2,3]):
  for n,color in [(2,'#00878a'),(8,'#b26126'),(32,'#77549b')]:
   rs=sorted([r for r in high if r['dimension']==dim and r['n']==n and r['curved']==1 and r['field']=='quintic' and r['method']=='fe_geometry'],key=lambda r:r['order'])
   ax.plot([r['order'] for r in rs],[max(1e-18,r['max_scaled_error']) for r in rs],'-o',markersize=3,color=color,label=f'h = 1/{n}')
  ax.axvspan(1,5,color='#00878a',alpha=.06);ax.axhline(1e-8,color='#777',linestyle=':',linewidth=1);ax.set_yscale('log');ax.set_xticks(range(1,11));ax.set_xlabel('Derivative order');ax.set_title('2D separable geometry' if dim==2 else '3D curved shear geometry');ax.grid(alpha=.2);ax.legend(fontsize=7)
 axes[0].set_ylabel('Max error / max(1, |reference|)');fig.savefig(figdir/'high-order.png',dpi=180);plt.close(fig)
 fig=plt.figure(figsize=(9,3.1),layout='constrained');ax=fig.add_subplot(121);az=fig.add_subplot(122,projection='3d')
 t=np.linspace(0,1,100)
 for s in np.linspace(0,1,9):
  b=1.2*s*(1-s)*t*(1-t);ax.plot(s+b,t+.7*b,color='#00878a',lw=.6);b=1.2*t*(1-t)*s*(1-s);ax.plot(t+b,s+.7*b,color='#00878a',lw=.6)
 ax.set_aspect('equal');ax.set_title('Exactly represented 2D coupled map');ax.set_xlabel('x');ax.set_ylabel('y')
 for d in range(3):
  for a in np.linspace(0,1,5):
   for b in [0,1]:
    coords=[np.full_like(t,a),np.full_like(t,b)];coords.insert(d,t);s,u,r=coords;x=s+.2*u*(1-u);y=u+.15*r*(1-r);az.plot(x,y,r,color='#00878a',lw=.7)
 az.set_title('Q2 curved 3D geometry');az.set_xlabel('x');az.set_ylabel('y');az.set_zlabel('z');az.set_box_aspect((1,1,1));fig.savefig(figdir/'meshes.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
