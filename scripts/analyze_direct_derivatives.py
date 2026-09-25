#!/usr/bin/env python3
"""Independent Decimal differentiation controls and matched-input summaries."""
import csv, json, math, statistics, time
from decimal import Decimal as D, localcontext
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/direct-derivatives'

def reference_derivative(mapping,nx,ny,step):
    # Independent physical inversion by Newton iteration at 65 decimal digits.
    # Central finite differences act on values, not on Taylor jets or derivatives.
    with localcontext() as context:
        context.prec=65
        a=D('1.2' if mapping=='coupled' else '.7')
        def geometry(s,t):
            if mapping=='coupled':
                b=a*s*(1-s)*t*(1-t)
                return s+b,t+D('.7')*b
            return s+a*s*(1-s),t+D('.7')*a*t*(1-t)
        s0=D('.12')+D('.20')*D('.23')+D('.06')*D('.31')
        t0=D('.18')-D('.03')*D('.23')+D('.22')*D('.31')
        x0,y0=geometry(s0,t0)
        def value(x,y):
            s,t=s0,t0
            for _ in range(14):
                gx,gy=geometry(s,t)
                if mapping=='coupled':
                    bs=a*(1-2*s)*t*(1-t);bt=a*s*(1-s)*(1-2*t)
                    j11,j12,j21,j22=1+bs,bt,D('.7')*bs,1+D('.7')*bt
                else:j11,j12,j21,j22=1+a*(1-2*s),D(0),D(0),1+D('.7')*a*(1-2*t)
                det=j11*j22-j12*j21;ex,ey=gx-x,gy-y
                s,t=s-(j22*ex-j12*ey)/det,t-(-j21*ex+j11*ey)/det
            ds,dt=s-D('.12'),t-D('.18');det=D('.20')*D('.22')+D('.06')*D('.03')
            r=(D('.22')*ds-D('.06')*dt)/det
            z=(D('.03')*ds+D('.20')*dt)/det
            return r**3*z+D('.4')*r*z*z+D('.7')*r*r-D('.2')*z
        stencils={0:[(0,D(1))],1:[(-1,D('-.5')),(1,D('.5'))],2:[(-1,D(1)),(0,D(-2)),(1,D(1))],3:[(-2,D('-.5')),(-1,D(1)),(1,D(-1)),(2,D('.5'))]}
        h=D(step)
        total=sum(wx*wy*value(x0+i*h,y0+j*h) for i,wx in stencils[nx] for j,wy in stencils[ny])
        return float(total/h**(nx+ny))

def main():
    start=time.monotonic()
    rows=list(csv.DictReader((OUT/'results.csv').open()))
    categorical={'family','basis','map','field','derivative','method'}
    for r in rows:
        for k in r:
            if k not in categorical:r[k]=float(r[k])
    key=lambda r:tuple(r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative'])
    bykey={key(r):r for r in rows if r['method']=='recursive'}
    comparisons=[]
    for r in rows:
        if r['method']!='direct':continue
        b=bykey[key(r)]
        result={k:r[k] for k in ['family','n','p','amplitude','basis','map','extra','field','derivative']}
        for metric in ['error_l2','face_rms','jump_rms']:
            result[metric+'_ratio']=r[metric]/b[metric] if b[metric]>1e-9 else math.nan
        comparisons.append(result)
    controls={k:max(r[k] for r in rows) for k in ['inverse_residual','representation_residual','first_derivative_residual','geometry_residual']}
    controls['constant_max_l2']=max(r['error_l2'] for r in rows if r['field']=='constant' and r['method']=='direct')
    controls['affine_quadratic_max_l2']=max(r['error_l2'] for r in rows if r['field']=='physical_quadratic' and r['method']=='direct')
    controls['represented_max_relative']=max(r['relative_l2'] for r in rows if r['field']=='represented_quadratic' and r['method']=='direct')
    lookup={key(r)+(r['method'],):r for r in rows}
    for name,selector in [('basis',lambda r:r['basis']!='gll'),('quadrature',lambda r:r['extra']>0)]:
        differences=[]
        for r in rows:
            if not selector(r):continue
            canonical=dict(r);canonical['basis']='gll';canonical['extra']=0
            b=lookup[key(canonical)+(r['method'],)]
            differences.append(abs(r['error_l2']-b['error_l2'])/max(1.,r['truth_l2']))
        controls[name+'_max_scaled_norm_change']=max(differences)
    previous=list(csv.DictReader((ROOT/'measurements/derivative-accuracy/results.csv').open()))
    diffs=[]
    for b in previous:
        if float(b['qextra'])!=0 or float(b['evalextra'])!=0:continue
        k=(b['family'],float(b['n']),float(b['p']),float(b['amplitude']),b['basis'],b['map'],0.,b['field'],b['derivative'])
        if k in bykey:diffs.append(abs(bykey[k]['error_l2']-float(b['error_l2']))/max(1.,float(b['truth_l2'])))
    controls['baseline_matches']=len(diffs);controls['baseline_max_scaled_norm_change']=max(diffs)
    probes=[]
    for m,nx,ny,v in csv.reader((OUT/'control-probes.csv').open()):
        coarse=reference_derivative(m,int(nx),int(ny),'1e-5');fine=reference_derivative(m,int(nx),int(ny),'1e-6')
        probes.append({'map':m,'nx':int(nx),'ny':int(ny),'jet':float(v),'decimal_fd_coarse':coarse,'decimal_fd_fine':fine,'scaled_difference':abs(float(v)-fine)/max(1.,abs(fine)),'reference_step_sensitivity':abs(coarse-fine)/max(1.,abs(fine))})
    controls['independent_fd_max_scaled_difference']=max(x['scaled_difference'] for x in probes)
    controls['independent_fd_max_step_sensitivity']=max(x['reference_step_sensitivity'] for x in probes)
    assert controls['independent_fd_max_scaled_difference']<1e-7,controls
    assert controls['baseline_max_scaled_norm_change']<1e-9,controls
    assert controls['represented_max_relative']<1e-6,controls
    representatives=[]
    for field,m,a in [('sine','coupled',1.2),('represented_quadratic','separable',.7)]:
        for family in ['quad','tri']:
            for n in [4,8,16]:
                for order in [1,2,3]:
                    for mode in ['recursive','direct']:
                        candidates=[r for r in rows if r['field']==field and r['family']==family and r['n']==n and r['p']==3 and r['basis']=='gll' and r['extra']==0 and r['amplitude']==a and len(r['derivative'])==order and r['method']==mode]
                        if candidates:representatives.append(max(candidates,key=lambda r:r['relative_l2']))
    aggregates=[]
    for m in ['coupled','separable']:
        selected=[r for r in comparisons if r['n']==8 and r['basis']=='gll' and r['extra']==0 and r['map']==m and r['amplitude']==(1.2 if m=='coupled' else .7) and r['field']!='constant']
        entry={'map':m,'count':len(selected)}
        for metric in ['error_l2','face_rms','jump_rms']:
            values=[r[metric+'_ratio'] for r in selected if math.isfinite(r[metric+'_ratio'])]
            entry[metric+'_median_ratio']=statistics.median(values);entry[metric+'_min_ratio']=min(values);entry[metric+'_max_ratio']=max(values);entry[metric+'_improved']=sum(v<1 for v in values)
        aggregates.append(entry)
    for name,data in [('comparisons',comparisons),('representatives',representatives),('independent-controls',probes)]:
        with (OUT/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    summary={'controls':controls,'aggregate':aggregates,'representatives':representatives,'row_count':len(rows),'independent_control_and_analysis_seconds':time.monotonic()-start}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    def find(field,family,n,mode):return next(r for r in representatives if r['field']==field and r['family']==family and r['n']==n and r['method']==mode and len(r['derivative'])==3)
    for family,color in [('quad','#007f86'),('tri','#c55420')]:
        label='Q3' if family=='quad' else 'P3'
        for mode,style,marker in [('recursive','-','o'),('direct','--','s')]:
            ns=[4,8,16];ys=[find('represented_quadratic',family,n,mode)['relative_l2'] for n in ns]
            axes[0].plot(ns,ys,style,color=color,marker=marker,label=label+' '+mode)
        ratio=[find('sine',family,n,'direct')['relative_l2']/find('sine',family,n,'recursive')['relative_l2'] for n in [4,8,16]]
        axes[1].plot([4,8,16],ratio,'-o',color=color,label=label)
    axes[0].set_yscale('log');axes[0].set_title('Exactly represented input');axes[0].set_ylabel('Worst third-derivative relative L2 error');axes[0].legend(fontsize=8)
    axes[1].axhline(1,color='#666666',linestyle=':',label='Equal errors');axes[1].set_title('Projected sine: direct / recursive');axes[1].set_ylabel('Ratio of worst third-derivative errors');axes[1].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('Subdivisions per direction, n');ax.set_xticks([4,8,16]);ax.grid(alpha=.2)
    (OUT/'figures').mkdir(exist_ok=True);fig.savefig(OUT/'figures/comparison.png',dpi=180);plt.close(fig)
    print(json.dumps({'controls':controls,'aggregate':aggregates},indent=2))
if __name__=='__main__':main()
