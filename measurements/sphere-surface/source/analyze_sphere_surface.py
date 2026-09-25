#!/usr/bin/env python3
"""Analyze the embedded-sphere study and render actual surface samples."""
import csv
from collections import Counter
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'measurements/sphere-surface'

def main():
    rows=list(csv.DictReader((OUT/'results.csv').open()))
    for row in rows:
        for k in row:
            if k!='family':row[k]=float(row[k])
    mainrows=[r for r in rows if r['p']==r['g'] and r['extra']==0 and r['constant']==0]
    finest=[]
    for family in ['tri','quad']:
        for p in [1,2,3]:
            sequence=sorted([r for r in mainrows if r['family']==family and r['p']==p],key=lambda r:r['levels'])
            coarse,fine=sequence[-2:];r=dict(fine);r['l2_rate']=math.log(coarse['l2']/fine['l2'],2);r['h1_rate']=math.log(coarse['h1_seminorm']/fine['h1_seminorm'],2);finest.append(r)
    sensitivity=[]
    for r in rows:
        if r['extra']:
            b=next(b for b in rows if b['family']==r['family'] and b['p']==r['p'] and b['g']==r['g'] and b['levels']==r['levels'] and b['extra']==0 and b['constant']==r['constant'])
            sensitivity.append(abs(r['l2']-b['l2'])/max(1.,b['l2']))
    controls={'constant_relative_l2_max':max(r['relative_l2'] for r in rows if r['constant']), 'quadrature_scaled_l2_change_max':max(sensitivity),'algebraic_relative_residual_max':max(r['relative_residual'] for r in rows),'mean_balance_defect_max':max(r['mean_balance_defect'] for r in rows)}
    assert controls['constant_relative_l2_max']<1e-10
    assert controls['quadrature_scaled_l2_change_max']<1e-9
    assert controls['algebraic_relative_residual_max']<1e-9
    assert all(r['l2_rate']>r['p']+.5 for r in finest),finest
    summary={'finest':finest,'controls':controls,'geometry_comparison':[r for r in rows if r['levels']==3 and r['p']==2 and r['extra']==0 and r['constant']==0]}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(9,3.1),layout='constrained')
    for ax,family in zip(axes,['tri','quad']):
        for p,color,marker in [(1,'#888888','o'),(2,'#007f86','s'),(3,'#c55420','^')]:
            values=sorted([r for r in mainrows if r['family']==family and r['p']==p],key=lambda r:r['levels'])
            ax.semilogy([r['levels'] for r in values],[r['relative_l2'] for r in values],'-',marker=marker,color=color,label=f'p=g={p}')
        ax.set_title('Triangles' if family=='tri' else 'Quadrilaterals');ax.set_xlabel('Uniform refinement level');ax.set_ylabel('Relative surface L2 error');ax.set_xticks([1,2,3,4]);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(OUT/'figures/convergence.png',dpi=180);plt.close(fig)
    # These vertices sample the actual approximate mesh Gamma_h. Do not radially
    # normalize them for rendering: that would hide geometry approximation.
    samples=np.genfromtxt(OUT/'visualization/tri-surface.csv',delimiter=',',names=True)
    vertices=np.column_stack([samples[k] for k in ['x','y','z']]).reshape(-1,3,3)
    solution=samples['u'].reshape(-1,3).mean(axis=1)
    error=samples['error'].reshape(-1,3).mean(axis=1)
    # Recover boundaries of original elements by cancelling shared rendering
    # edges within each element. Subdivision edges are not mesh edges.
    element=samples['element'][::3].astype(int)
    edges=[]
    for e in np.unique(element):
        counts=Counter()
        for tri in vertices[element==e]:
            for i,j in [(0,1),(1,2),(2,0)]:
                key=tuple(sorted((tuple(np.round(tri[i],12)),tuple(np.round(tri[j],12)))))
                counts[key]+=1
        edges.extend(key for key,count in counts.items() if count==1)
    view=np.array([math.cos(math.radians(23))*math.cos(math.radians(35)),math.cos(math.radians(23))*math.sin(math.radians(35)),math.sin(math.radians(23))])
    edges=[edge for edge in edges if np.dot(np.mean(edge,axis=0),view)>0.12]
    fig=plt.figure(figsize=(9,3.7));norms=[Normalize(-.5,.5),Normalize(-np.max(np.abs(samples['error'])),np.max(np.abs(samples['error'])))];cmaps=[plt.get_cmap('coolwarm'),plt.get_cmap('RdBu_r')]
    for i,(values,title) in enumerate([(solution,'Computed solution on the curved mesh'),(error,'Signed error against radial exact extension')]):
        ax=fig.add_subplot(1,2,i+1,projection='3d',computed_zorder=False);surface=Poly3DCollection(vertices,facecolors=cmaps[i](norms[i](values)),edgecolors='none',linewidths=0,antialiased=False,rasterized=True);surface.set_zorder(1);ax.add_collection3d(surface)
        if i==0:ax.add_collection3d(Line3DCollection(edges,colors='#263840',linewidths=.45,alpha=.8,zorder=2))
        ax.set(xlim=(-1.05,1.05),ylim=(-1.05,1.05),zlim=(-1.05,1.05));ax.set_box_aspect((1,1,1));ax.view_init(elev=23,azim=35);ax.set_axis_off();ax.set_title(title,fontsize=10,pad=1)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norms[i],cmap=cmaps[i]),ax=ax,fraction=.035,pad=.01,shrink=.75);cb.ax.tick_params(labelsize=8);cb.set_label('u_h' if i==0 else 'u_h - u_exact',fontsize=9)
        if i:cb.formatter.set_powerlimits((0,0));cb.update_ticks()
    fig.subplots_adjust(left=.005,right=.94,top=.91,bottom=.04,wspace=.09);fig.savefig(OUT/'figures/sphere.png',dpi=200);plt.close(fig)
    fig=plt.figure(figsize=(9,1.05));fig.text(.02,.65,r'$G=J^{T}J,\qquad \nabla_{\Gamma_h}v=JG^{-1}\widehat{\nabla}\widehat{v},\qquad dS=\sqrt{\det G}\,d\xi$',fontsize=13)
    fig.text(.02,.13,r'$\int_{\Gamma_h}(\nabla_{\Gamma_h}u_h\cdot\nabla_{\Gamma_h}v_h+u_hv_h)\,dS=\int_{\Gamma_h}f^{e}v_h\,dS$',fontsize=13)
    fig.savefig(OUT/'figures/surface-math.png',dpi=200);plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
