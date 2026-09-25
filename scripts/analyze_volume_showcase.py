"""Generate diagnostic figures from retained 3D polynomial samples."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; D=ROOT/'measurements/volume-showcase'; F=D/'figures';F.mkdir(exist_ok=True)
rows=list(csv.DictReader((D/'results.csv').open()))
for r in rows:
    for k in r:
        if k!='method':r[k]=float(r[k])
(D/'summary.json').write_text(json.dumps({'rows':rows},indent=2))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def data(method,state):return np.genfromtxt(D/'visualization'/f'n10-{method}-p1-dt0.003-{state}.csv',delimiter=',',names=True)
dg=data('dg','final');oe=data('ofdg','final');ini=data('ofdg','initial')
fig=plt.figure(figsize=(10,3.5))
for i,(a,label) in enumerate([(ini,'Initial OFDG field'),(dg,'Final DG'),(oe,'Final OFDG')]):
    ax=fig.add_subplot(1,3,i+1,projection='3d');mask=a['u']>.5
    ax.scatter(a['x'][mask],a['y'][mask],a['z'][mask],s=5,alpha=.6,color=['#63717d','#3366a8','#df7327'][i])
    ax.set(xlim=(0,1),ylim=(0,1),zlim=(0,1),xlabel='x',ylabel='y',zlabel='z',title=label);ax.set_box_aspect((1,1,1));ax.view_init(24,-57)
fig.tight_layout();fig.savefig(F/'volume.png',dpi=190);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(10,3.2),layout='constrained')
lo=min(dg['u'].min(),oe['u'].min());hi=max(dg['u'].max(),oe['u'].max())
for ax,a,title in zip(axes,[dg,oe,oe],['DG: z = 0.45','OFDG: z = 0.45','Exact: z = 0.45']):
    m=np.isclose(a['z'],.45);b=a[m];ix=np.rint(b['x']*50-.5).astype(int);iy=np.rint(b['y']*50-.5).astype(int);v=np.empty((50,50));v[iy,ix]=b['exact' if title.startswith('Exact') else 'u'];assert len(set(zip(ix,iy)))==2500
    im=ax.imshow(v,origin='lower',extent=(0,1,0,1),vmin=lo,vmax=hi,cmap='viridis');ax.set(xlabel='x',ylabel='y',title=title,xlim=(.1,.8),ylim=(.1,.8))
fig.colorbar(im,ax=axes,shrink=.8,label='u');fig.savefig(F/'slice.png',dpi=190);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,2.8),layout='constrained')
for ax in axes:
    for a,label,color in [(dg,'DG','#3366a8'),(oe,'OFDG','#df7327')]:
        m=np.isclose(a['y'],.43)&np.isclose(a['z'],.45);b=a[m];b=b[np.argsort(b['x'])]
        # Break between elements instead of connecting discontinuous DG traces.
        for e in np.unique(b['element']):
            c=b[b['element']==e];ax.plot(c['x'],c['u'],'.-',color=color,label=label if e==np.unique(b['element'])[0] else None,ms=3)
    radius=np.sqrt(.04-(.43-.44)**2-(.45-.46)**2);left=.45-radius;right=.45+radius
    ax.plot([0,left,left,right,right,1],[0,0,1,1,0,0],'k--',lw=1,label='Exact')
    ax.set(xlabel='x',ylabel='u');ax.grid(alpha=.2)
axes[0].set(title='Physical line: y=0.43, z=0.45',xlim=(.1,.8));axes[0].legend(ncol=3,fontsize=9)
axes[1].set(title='Enlarged leading edge',xlim=(.57,.74),ylim=(-.35,1.4))
fig.savefig(F/'cut.png',dpi=190);plt.close(fig)
