#!/usr/bin/env python3
"""Analyze the isolated first-order surface-OFDG experiment."""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/sphere-ofdg'

def main():
    rows=list(csv.DictReader((OUT/'results.csv').open()))
    for r in rows:
        for k in r:
            if k not in ['family','method']:r[k]=float(r[k])
    base=[r for r in rows if r['cfl']==.12 and r['quadrature']==20]
    def select(family,level,method,problem,cfl=.12,q=20):return next(r for r in rows if r['family']==family and r['level']==level and r['method']==method and r['problem']==problem and r['cfl']==cfl and r['quadrature']==q)
    finest=[]
    for family in ['tri','quad']:
        for method in ['dg','ofdg']:
            r=dict(select(family,3,method,0));r['l2_rate']=math.log(select(family,2,method,0)['l2']/r['l2'],2);finest.append(r)
    controls={key:max(r[key] for r in rows) for key in ['mass_drift','filter_mean_defect','constant_operator_residual','flux_mismatch']}
    controls['constant_error_l2']=max(r['l2'] for r in rows if r['problem']==2)
    controls['area_error']=max(abs(r['area']-4*math.pi) for r in rows)
    controls['time_sensitivity']={m:abs(select('tri',3,m,0,cfl=.06)['l2']/select('tri',3,m,0)['l2']-1) for m in ['dg','ofdg']}
    controls['smooth_quadrature_sensitivity']={m:abs(select('tri',2,m,0,q=26)['l2']/select('tri',2,m,0)['l2']-1) for m in ['dg','ofdg']}
    controls['cap_quadrature_sensitivity']={m:abs(select('tri',3,m,1,q=30)['l2']/select('tri',3,m,1)['l2']-1) for m in ['dg','ofdg']}
    assert all(math.isfinite(r['l2']) and r['l2']>=0 for r in rows)
    assert controls['constant_error_l2']<1e-10 and controls['mass_drift']<1e-9
    assert controls['filter_mean_defect']<1e-12 and controls['constant_operator_residual']<1e-8
    summary={'finest_smooth':finest,'finest_cap':[select(f,3,m,1) for f in ['tri','quad'] for m in ['dg','ofdg']], 'controls':controls}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(9,3),layout='constrained')
    for ax,family in zip(axes,['tri','quad']):
        for method,color,marker in [('dg','#6c7886','o'),('ofdg','#007f86','s')]:
            rr=[select(family,l,method,0) for l in [1,2,3]];ax.loglog([r['dofs'] for r in rr],[r['relative_l2'] for r in rr],'-',color=color,marker=marker,label=method.upper())
        ax.set_title('P1 triangles' if family=='tri' else 'Q1 quadrilaterals');ax.set_xlabel('Scalar unknowns');ax.set_ylabel('Relative L2 error');ax.set_xticks([96,384,1536],labels=['96','384','1536']);ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.savefig(OUT/'figures/smooth.png',dpi=180);plt.close(fig)
    datasets=[np.genfromtxt(OUT/'visualization'/name,delimiter=',',names=True) for name in ['tri-dg-problem1-initial.csv','tri-dg-problem1-final.csv','tri-ofdg-problem1-final.csv']]
    bound=max(.2,max(np.max(d['u'])-1 for d in datasets[1:]),max(-np.min(d['u']) for d in datasets[1:]))
    norm=Normalize(-bound,1+bound);cmap=plt.get_cmap('viridis')
    fig=plt.figure(figsize=(9,3.1))
    for i,(d,title) in enumerate(zip(datasets,['Exact initial cap','DG after quarter rotation','Surface OFDG after quarter rotation'])):
        ax=fig.add_subplot(1,3,i+1,projection='3d');v=np.column_stack([d[k] for k in ['x','y','z']]).reshape(-1,3,3);values=d['exact' if i==0 else 'u'].reshape(-1,3).mean(axis=1)
        ax.add_collection3d(Poly3DCollection(v,facecolors=cmap(norm(values)),edgecolors='none',linewidths=0,antialiased=False,rasterized=True));ax.set(xlim=(-1.02,1.02),ylim=(-1.02,1.02),zlim=(-1.02,1.02));ax.set_box_aspect((1,1,1));ax.view_init(20,45);ax.set_axis_off();ax.set_title(title,fontsize=9,pad=-3)
    fig.subplots_adjust(left=0,right=.88,top=.87,bottom=.01,wspace=-.04);cax=fig.add_axes([.91,.18,.014,.60]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.ax.tick_params(labelsize=8);cb.set_label('u (shared scale)',fontsize=8);fig.savefig(OUT/'figures/cap-sphere.png',dpi=190);plt.close(fig)
    # Equator is an element boundary. Use northern one-sided traces, and keep
    # separate elements as separate lines rather than smoothing DG jumps.
    fig,axes=plt.subplots(1,2,figsize=(9,2.9),layout='constrained')
    for ax,limits in zip(axes,[(-20,200),(10,55)]):
        theta=np.array([-20,30,30,150,150,200]);exact=np.array([0,0,1,1,0,0]);ax.plot(theta,exact,'k--',lw=1.2,label='Exact')
        for data,method,color,marker in [(datasets[1],'DG','#6c7886','o'),(datasets[2],'OFDG','#007f86','s')]:
            first=True
            for e in np.unique(data['element']):
                block=data[data['element']==e]
                if np.mean(block['z'])<=1e-12:continue
                pts=block[np.abs(block['z'])<1e-12]
                if len(pts)<2:continue
                theta=np.degrees(np.arctan2(pts['y'],pts['x']));anchor=np.degrees(np.arctan2(np.mean(pts['y']),np.mean(pts['x'])));theta=anchor+(theta-anchor+180)%360-180
                if anchor<-20:theta+=360
                values=np.unique(np.column_stack([theta,pts['u']]),axis=0);values=values[np.argsort(values[:,0])]
                if values[-1,0]<limits[0] or values[0,0]>limits[1]:continue
                ax.plot(values[:,0],values[:,1],color=color,marker=marker,markersize=2.3,linewidth=1.1,label=method if first else None);first=False
        ax.set_xlim(*limits);ax.set_ylim(-.23,1.23);ax.set_xlabel('Equatorial longitude (degrees)');ax.set_ylabel('u');ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes[0].set_title('Final equatorial trace');axes[1].set_title('Cap edge: enlarged view');fig.savefig(OUT/'figures/cap-cut.png',dpi=190);plt.close(fig)
    fig=plt.figure(figsize=(9,1.8))
    equations=[r'$g_{K,a}=\Pi_K(\nabla_\Gamma u_h)_a,\qquad A=\max|u_h-\bar u|,\quad h_K=(|K|/|\widehat K|)^{1/2}$',r'$J_{F,0}=|F|^{-1}\int_F|[u_h]|\,ds,\quad J_{F,1}=|F|^{-1}\int_F\sum_{a=1}^{3}|[g_a]|\,ds$',r'$\lambda_K=A^{-1}\sum_{F\subset\partial K}\beta_F\left(\frac{J_{F,0}}{2h_K}+\frac{3J_{F,1}}{2}\right),\quad u_h^+=\bar u_K+e^{-\Delta t\lambda_K}(u_h^- -\bar u_K)$']
    for i,eq in enumerate(equations):fig.text(.015,.77-i*.32,eq,fontsize=12)
    fig.savefig(OUT/'figures/formula.png',dpi=200);plt.close(fig)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
