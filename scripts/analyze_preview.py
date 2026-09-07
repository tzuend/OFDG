#!/usr/bin/env python3
"""Create preview-only figures and tables from checkpointed measured results."""
import argparse,csv,json,math,os,subprocess,sys
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/ofdg-preview-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from matplotlib.lines import Line2D
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
METHODS=['dg','ofdg','ofdg-kxrcf','oedg']
LABELS={'dg':'DG','ofdg':'OFDG','ofdg-kxrcf':'OFDG-KXRCF','oedg':'OEDG'}
COLORS=dict(zip(METHODS,['#697586','#e07a27','#007f86','#8656a2']))
STYLES={'dg':'-', 'ofdg':'--', 'ofdg-kxrcf':':', 'oedg':'-.'}
MARKERS={'dg':'s', 'ofdg':'D', 'ofdg-kxrcf':'o', 'oedg':'^'}
plt.rcParams.update({'font.size':9.5, 'axes.titlesize':10, 'axes.labelsize':9.5,
    'xtick.labelsize':9, 'ytick.labelsize':9, 'axes.spines.top':False,
    'axes.spines.right':False, 'axes.grid':True, 'grid.alpha':.16,
    'lines.linewidth':1.35, 'legend.fontsize':9, 'pdf.fonttype':42})


def method_style(method):
    return dict(color=COLORS[method], linestyle=STYLES[method],
                linewidth=1.6 if method=='dg' else 1.25)


def method_legend(fig, methods=METHODS, reference=False):
    handles=[Line2D([],[],**method_style(m),marker=MARKERS[m],markersize=4,
                    markerfacecolor='white',label=LABELS[m]) for m in methods]
    if reference:handles.append(Line2D([],[],color='black',linestyle=(0,(3,2)),linewidth=1,label='Reference'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.015),
               ncol=len(handles),frameon=False,handlelength=2.7,columnspacing=1.2)


def lagrange_weights(nodes, positions):
    nodes=np.asarray(nodes); positions=np.atleast_1d(positions)
    basis=np.ones((len(positions),len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i!=j:basis[:,i]*=(positions-nodes[j])/(nodes[i]-nodes[j])
    return basis



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'measurements/study/preview')
    parser.add_argument('--pdf-python',default=str(Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3'))
    parser.add_argument('--output',type=Path,help='Separate analysis destination; raw inputs remain unchanged')
    parser.add_argument('--report',action='store_true',help='Narrower thesis figures; do not rebuild the supervisor brief')
    args=parser.parse_args(); directory=args.input.resolve()
    if args.output and not args.report:
        parser.error('--output currently requires --report to keep the supervisor brief separate')
    output=(args.output or (ROOT/'report/generated/preview' if args.report else directory)).resolve()
    if args.report and output==directory:
        parser.error('Report output must be separate from raw preview inputs')
    output.mkdir(parents=True,exist_ok=True)
    state=json.loads((directory/'results.json').read_text())
    figures=output/'figures';figures.mkdir(exist_ok=True)
    runs=state['runs']; success=[r for r in runs if r['status']=='ok']
    def select(case,method=None,n=None,p=None):
        selected=[r for r in success if r['case']==case and (method is None or r['method']==method) and (n is None or r['n']==n) and (p is None or r['order']==p)]
        # If a sensitivity check triggered a rerun, use the smaller CFL consistently.
        corrected=case=='smooth_advection' and p==3 and state.get('temporal_sensitivity',{}).get(method,0)>.05
        cfl=.075 if corrected else None
        return sorted([r for r in selected if r['cfl']==(cfl or (.3 if r['program']=='euler' else .15))],key=lambda r:r['n'])
    cache={}
    def samples(run,kind='polynomial'):
        key=(run['id'],kind)
        if key not in cache:
            data=[]
            for path in directory.joinpath('profiles').glob(run['id']+'.rank*.csv'):
                for row in csv.DictReader(path.open()):
                    if row['sample']==kind:
                        data.append({k:float(v) for k,v in row.items() if k!='sample'})
            cache[key]=sorted(data,key=lambda r:(r['x'],r['y'],r['rank'],r['element']))
        return cache[key]
    def savefig(fig,name):
        if args.report:
            fig.set_size_inches(5.1, fig.get_figheight()+.35)
            for legend in list(fig.legends):
                legend.remove()
            if name not in ['riemann', 'nonlinear_curved_accuracy']:
                methods=METHODS if name in ['scalar_shocks','euler_profiles'] else METHODS[:3]
                handles=[Line2D([],[],**method_style(m),marker=MARKERS[m],markersize=3,markerfacecolor='white',label=LABELS[m]) for m in methods]
                if name in ['scalar_shocks','euler_profiles']:
                    handles.append(Line2D([],[],color='black',ls='--',label='Reference'))
                fig.legend(handles=handles,loc='upper center',ncol=3,frameon=False,fontsize=8.5)
            if name!='riemann':
                fig.tight_layout(rect=(0,0,1,.89 if name in ['accuracy','riemann_slice'] else .92))
        fig.savefig(figures/(name+'.png'),dpi=280,bbox_inches='tight')
        fig.savefig(figures/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    convergence=[]
    for case in ['smooth_advection','smooth_burgers','straight_advection_2d','curved_advection_2d']:
        for method in METHODS[:3]:
            for p in ([1,2,3] if case=='smooth_advection' else [2]):
                series=select(case,method,p=p)
                for i,r in enumerate(series):
                    error=float(r['metrics']['l2_error'])
                    order=math.log(float(series[i-1]['metrics']['l2_error'])/error)/math.log(r['n']/series[i-1]['n']) if i and error>0 else None
                    convergence.append(dict(case=case,method=method,p=p,n=r['n'],l2=error,observed_order=order,cfl=r['cfl']))
    with (output/'convergence.csv').open('w') as f:
        w=csv.DictWriter(f,list(convergence[0]));w.writeheader();w.writerows(convergence)
    fig,axes=plt.subplots(1,3,figsize=(7.4,2.55))
    for p,ax in zip([1,2,3],axes):
        for method in METHODS[:3]:
            rows=select('smooth_advection',method,p=p)
            ax.loglog([r['n'] for r in rows],[float(r['metrics']['l2_error']) for r in rows],
                      **method_style(method),marker=MARKERS[method],
                      ms=6 if method=='dg' else 3.8,mfc='white',mew=1.1)
        ax.set(title=f'P{p} advection',xlabel='Cells',ylabel='L2 error')
        ax.set_xticks([32,64,128,256],labels=['32','64','128','256'])
        ax.xaxis.set_minor_formatter(NullFormatter())
    method_legend(fig,METHODS[:3]);fig.tight_layout(rect=(0,0,1,.91));savefig(fig,'accuracy')
    fig,axes=plt.subplots(1,3,figsize=(7.4,2.45))
    for ax,case,title in zip(axes,['smooth_burgers','straight_advection_2d','curved_advection_2d'],
                            ['P2 Burgers, T=0.6','P2 straight, T=1','P2 curved, T=1']):
        for method in METHODS[:3]:
            rows=select(case,method,p=2)
            ax.loglog([r['n'] for r in rows],[float(r['metrics']['l2_error']) for r in rows],
                      **method_style(method),marker=MARKERS[method],
                      ms=6 if method=='dg' else 3.8,mfc='white',mew=1.1)
        ticks=[32,64,128] if case=='smooth_burgers' else [8,16,32]
        ax.set(title=title,xlabel='Cells' if case=='smooth_burgers' else 'Cells / direction',ylabel='L2 error')
        ax.set_xticks(ticks,labels=list(map(str,ticks)));ax.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout();savefig(fig,'nonlinear_curved_accuracy')
    def plot_profile(ax,case,n,p=2,field='value',methods=METHODS,window=None):
        for method in methods:
            rows=select(case,method,n,p)
            if not rows:continue
            r=rows[0];data=samples(r)
            # Break curves between cells: DG jumps must remain visible.
            groups={}
            for s in data:groups.setdefault((s['rank'],s['element']),[]).append(s)
            x=[];y=[]
            for points in groups.values():
                points.sort(key=lambda s:s['x']);x.extend([s['x'] for s in points]+[np.nan]);y.extend([s[field] for s in points]+[np.nan])
            ax.plot(x,y,**method_style(method),label=LABELS[method])
            if window:
                visible=[point for point in data if window[0]<=point['x']<=window[1]]
                spacing=max(1,len(visible)//8)
                offset=METHODS.index(method)*max(1,spacing//4)
                marked=visible[offset::spacing]
                ax.plot([v['x'] for v in marked],[v[field] for v in marked],linestyle='none',
                        marker=MARKERS[method],ms=3.8,mfc='white',mew=.9,color=COLORS[method])
                ax.set_xlim(window)
        ax.set(xlabel='x',ylabel=field)
    def filtered_difference(ax,case,n,window,field='value'):
        base=samples(select(case,'ofdg',n,2)[0])
        maxima={}
        for method in ['ofdg-kxrcf','oedg']:
            other=samples(select(case,method,n,2)[0])
            assert len(base)==len(other)
            xx=[];dd=[];last=None
            for a,b in zip(base,other):
                assert abs(a['x']-b['x'])<1e-12
                if not window[0]<=a['x']<=window[1]:continue
                key=(a['rank'],a['element'])
                if last is not None and key!=last:xx.append(np.nan);dd.append(np.nan)
                xx.append(a['x']);dd.append(b[field]-a[field]);last=key
            maxima[method]=float(np.nanmax(np.abs(dd)))
            ax.plot(xx,dd,**method_style(method),marker=MARKERS[method],ms=3,
                    mfc='white',markevery=max(1,len(xx)//7))
        ax.axhline(0,color=COLORS['ofdg'],lw=.8,ls='--')
        ax.set(xlim=window,xlabel='x',ylabel='u - u(OFDG)')
        ax.yaxis.set_major_locator(MaxNLocator(4))
        return maxima
    scalar_differences={}
    fig,axes=plt.subplots(3,2,figsize=(7.4,5.65))
    x=np.linspace(0,1,12000);foot=(x-1.1)%1
    exact=np.where((foot>=.3)&(foot<=.8),np.sin(2*np.pi*foot),np.cos(2*np.pi*foot)-.5)
    exact[np.abs(np.diff(exact,prepend=exact[0]))>.1]=np.nan
    center=math.pi+.5*2.2
    windows=[(.375,.425),(center-.13,center+.13)]
    for col,case in enumerate(['piecewise_advection','shock_burgers']):
        window=windows[col]
        for row in [0,1]:
            plot_profile(axes[row,col],case,256,window=window if row else None)
            if col==0:axes[row,col].plot(x,exact,color='black',ls=(0,(3,2)),lw=.8)
        axes[0,col].axvspan(*window,color='#c8e7e8',alpha=.35)
        axes[1,col].set_title('Jump detail' if col==0 else 'Shock detail')
        scalar_differences[case]=filtered_difference(axes[2,col],case,256,window)
        axes[2,col].set_title('Filtered-method difference')
    axes[0,0].set_title('Transport, P2 / 256 cells')
    axes[0,1].set_title('Burgers, P2 / 256 cells')
    method_legend(fig,reference=True);fig.tight_layout(rect=(0,0,1,.94));savefig(fig,'scalar_shocks')
    extrema=[]
    for case in ['piecewise_advection','shock_burgers']:
        lower,upper=(-1.,math.sin(.6*math.pi)) if case=='piecewise_advection' else (-.5,1.5)
        for method in METHODS:
            rows=select(case,method,256,2)
            if not rows:continue
            values=[s['value'] for s in samples(rows[0])]
            extrema.append(dict(case=case,method=method,minimum=min(values),maximum=max(values),undershoot=max(0,lower-min(values)),overshoot=max(0,max(values)-upper)))
    references={};sensitivity={}
    for case in ['lax','shu_osher']:
        for n in [2048,4096]:
            path=directory/'references'/f'{case}_{n}.csv'
            spec=next((r for r in state['references'] if r['name']==f'{case}_{n}'),{})
            if path.exists() and spec.get('status')=='ok':
                references[(case,n)]=np.genfromtxt(path,delimiter=',',names=True)
        if (case,2048) in references and (case,4096) in references:
            coarse=references[(case,2048)]['density'];fine=references[(case,4096)]['density'].reshape(-1,2).mean(axis=1)
            sensitivity[case]=dict(mean_absolute_density_difference=float(np.mean(abs(coarse-fine))),max_density_difference=float(np.max(abs(coarse-fine))))
    shock_centers={}
    for case in ['lax','shu_osher']:
        ref=references.get((case,4096))
        if ref is not None:
            mids=(ref['x'][1:]+ref['x'][:-1])/2
            gradient=np.abs(np.diff(ref['density']))
            if case=='shu_osher':gradient=np.where((mids>2.2)&(mids<2.8),gradient,0)
            shock_centers[case]=float(mids[np.argmax(gradient)])
        else:shock_centers[case]=3. if case=='lax' else 2.4
    wave_window=(.9,1.65); wave_ranges={}
    fig,axes=plt.subplots(3,2,figsize=(7.4,5.55))
    for col,(case,n) in enumerate([('lax',128),('shu_osher',400)]):
        center=shock_centers[case];half=.24 if case=='lax' else .13
        window=(center-half,center+half)
        plot_profile(axes[0,col],case,n,field='density')
        plot_profile(axes[1,col],case,n,field='density',window=window)
        axes[0,col].axvspan(*window,color='#c8e7e8',alpha=.35)
        axes[1,col].set_title('Shock-front detail')
        if (case,4096) in references:
            ref=references[(case,4096)]
            for row in [0,1]:axes[row,col].plot(ref['x'],ref['density'],color='black',ls=(0,(3,2)),lw=.85)
    axes[0,0].set_title('Lax, P2 / 128 cells')
    axes[0,1].set(title='Shu-Osher: cropped overview',xlim=(-1,3))
    axes[0,1].axvspan(*wave_window,color='#fae5cf',alpha=.45)
    plot_profile(axes[2,1],'shu_osher',400,field='density',window=wave_window)
    axes[2,1].set(title='Wave detail: x = 0.90 to 1.65',ylim=(3.0,4.85))
    if ('shu_osher',4096) in references:
        ref=references[('shu_osher',4096)]
        axes[2,1].plot(ref['x'],ref['density'],color='black',ls=(0,(3,2)),lw=.85)
        vals=ref['density'][(ref['x']>=wave_window[0])&(ref['x']<=wave_window[1])]
        wave_ranges['reference']=float(vals.max()-vals.min())
    for method in METHODS:
        vals=[d['density'] for d in samples(select('shu_osher',method,400,2)[0]) if wave_window[0]<=d['x']<=wave_window[1]]
        wave_ranges[method]=max(vals)-min(vals)
    # Equal-DOF comparison retains physical cell means, with explicit point markers.
    ax=axes[2,0];window=(shock_centers['lax']-.24,shock_centers['lax']+.24)
    for p,n,style,marker,color in [(1,192,'--','s','#b5611b'),(2,128,'-','o','#007f86'),(3,96,':','^','#8656a2')]:
        data=samples(select('lax','ofdg-kxrcf',n,p)[0],'cell_average')
        ax.plot([d['x'] for d in data],[d['density'] for d in data],style,color=color,marker=marker,ms=3,mfc='white',label=f'P{p}')
    ref=references.get(('lax',4096))
    if ref is not None:ax.plot(ref['x'],ref['density'],color='black',ls=(0,(3,2)),lw=.85)
    ax.set(xlim=window,title='Lax: equal-DOF cell means',xlabel='x',ylabel='Mean density')
    ax.legend(loc='upper right',ncol=3,fontsize=8.5,columnspacing=.6,handlelength=1.3,framealpha=.95)
    method_legend(fig,reference=True);fig.tight_layout(rect=(0,0,1,.94));savefig(fig,'euler_profiles')
    # Prefer the larger common successful resolution; never replace a failed method silently.
    common=[n for n in [32,64] if all(select('riemann_2d',m,n,2) for m in METHODS[:3])]
    n=max(common) if common else 32
    map_data={}
    for method in METHODS[:3]:
        data=samples(select('riemann_2d',method,n,2)[0],'cell_average')
        ordered=sorted(data,key=lambda d:(d['y'],d['x']))
        map_data[method]=np.array([d['density'] for d in ordered]).reshape(n,n)
    map_min=min(z.min() for z in map_data.values());map_max=max(z.max() for z in map_data.values())
    fig,axes=plt.subplots(1,3,figsize=(7.4,2.45),layout='constrained')
    for ax,method in zip(axes,METHODS[:3]):
        im=ax.imshow(map_data[method],origin='lower',extent=(0,1,0,1),vmin=map_min,vmax=map_max,
                     cmap='viridis',interpolation='nearest')
        ax.set(title=LABELS[method],xlabel='x',ylabel='y')
        ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    fig.colorbar(im,ax=axes,shrink=.8,label='Mean density')
    savefig(fig,'riemann')

    def density_cut(run,physical_y=.74):
        # Cartesian P2 tensor nodes determine the polynomial exactly. Interpolate
        # each cell independently; do not smooth or average across DG interfaces.
        centers=samples(run,'cell_average');nodes=samples(run)
        selected={(v['rank'],v['element']):v for v in centers
                  if abs(v['y']-physical_y)<.5/run['n']-1e-12}
        groups={key:[] for key in selected}
        for point in nodes:
            key=(point['rank'],point['element'])
            if key in groups:groups[key].append(point)
        x_all=[];rho_all=[]
        for key in sorted(selected,key=lambda k:selected[k]['x']):
            points=groups[key]
            # Bilinear coordinate evaluation can differ by one ulp along a tensor
            # row. Cluster only this roundoff; retain averaged physical locations.
            def coordinate_nodes(axis):
                clusters={}
                for point in points:clusters.setdefault(round(point[axis],12),[]).append(point[axis])
                return sorted(float(np.mean(values)) for values in clusters.values())
            xs=coordinate_nodes('x');ys=coordinate_nodes('y')
            assert len(xs)==len(ys)==run['order']+1
            data=np.array([[next(v['density'] for v in points if abs(v['x']-x)<1e-12 and abs(v['y']-y)<1e-12) for x in xs] for y in ys])
            assert np.max(abs(lagrange_weights(ys,ys)@data@lagrange_weights(xs,xs).T-data))<1e-12
            xx=np.linspace(selected[key]['x']-.5/run['n'],selected[key]['x']+.5/run['n'],21)
            rho=(lagrange_weights(ys,[physical_y])@data@lagrange_weights(xs,xx).T).ravel()
            x_all.extend([*xx,np.nan]);rho_all.extend([*rho,np.nan])
        assert len(selected)==run['n']
        return np.array(x_all),np.array(rho_all)
    # Independent interpolation check against an exactly represented tensor polynomial.
    test_nodes=np.array([.1,.4,.8]);test_grid=test_nodes[:,None]**2+test_nodes[None,:]**2+test_nodes[:,None]*test_nodes[None,:]
    assert abs((lagrange_weights(test_nodes,[.74])@test_grid@lagrange_weights(test_nodes,[.37]).T)[0,0]-(.74**2+.37**2+.74*.37))<1e-13
    cut_differences={}
    fig,axes=plt.subplots(2,2,figsize=(7.4,3.75))
    for col,size in enumerate([32,64]):
        cuts={m:density_cut(select('riemann_2d',m,size,2)[0]) for m in METHODS[:3]}
        cut_differences[size]={m:float(np.nanmax(abs(cuts[m][1]-cuts['ofdg'][1]))) for m in ['dg','ofdg-kxrcf']}
        for method,(xx,rho) in cuts.items():
            for row in [0,1]:
                axes[row,col].plot(xx,rho,**method_style(method))
                if row:
                    visible=np.where((xx>=.85)&(xx<=.95))[0]
                    marked=visible[METHODS.index(method)*2::max(1,len(visible)//8)]
                    axes[row,col].plot(xx[marked],rho[marked],linestyle='none',marker=MARKERS[method],ms=3.5,mfc='white',color=COLORS[method])
        axes[0,col].set_title(f'{size} x {size} cells, y = 0.740')
        axes[0,col].axvspan(.85,.95,color='#c8e7e8',alpha=.35)
        axes[1,col].set(title='Shock detail',xlim=(.85,.95))
    for ax in axes.flat:ax.set(xlabel='x',ylabel='Density')
    method_legend(fig,METHODS[:3]);fig.tight_layout(rect=(0,0,1,.92));savefig(fig,'riemann_slice')
    summary=dict(scalar_differences=scalar_differences,wave_ranges=wave_ranges,cut_differences=cut_differences,
                 riemann_map_range=[float(map_min),float(map_max)],convergence=convergence,extrema=extrema,reference_sensitivity=sensitivity,riemann_resolution=n,
                 failures=[{k:r[k] for k in ['id','status']} for r in runs if r['status']!='ok'],
                 periodic_conservation_max=max(float(r['metrics'].get('conservation_drift',0)) for r in success if r['program']!='euler'),
                 elapsed_minutes=state['elapsed_seconds']/60,peak_mib=max([r.get('peak_mib',0) for r in runs+state['references']]),
                 success_count=len(success),run_count=len(runs),temporal_sensitivity=state.get('temporal_sensitivity',{}))
    (output/'analysis.json').write_text(json.dumps(summary,indent=2)+'\n')
    pdf_python=args.pdf_python if Path(args.pdf_python).exists() else sys.executable
    if not args.report:
        subprocess.run([pdf_python,str(ROOT/'scripts/build_preview_pdf.py'),'--input',str(directory)],check=True)

if __name__=='__main__':main()
