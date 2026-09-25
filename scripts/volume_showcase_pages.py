"""Two-page showcase of the unchanged volume OFDG method."""
import json
from reportlab.platypus import Paragraph,Image,PageBreak,Spacer,Table,TableStyle
from reportlab.lib import colors

def append_volume_pages(story,directory,styles):
    rows=json.loads((directory/'summary.json').read_text())['rows'];run=json.loads((directory/'run.json').read_text())
    def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
    def title(n,s):story.append(PageBreak());text(f'3D VOLUME SHOWCASE / {n}','SmallPreview');text(s,'TitlePreview')
    def image(name,height):story.append(Image(str(directory/'figures'/name),width=504,height=height));story.append(Spacer(1,5))
    def table(data,widths):
        t=Table([[Paragraph(str(x),styles['TablePreview']) for x in r] for r in data],colWidths=widths)
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,7)])
    title(18,'A genuine three-dimensional volume test')
    text('<b>The existing production OFDG implementation runs this example without modification.</b> The domain is the volume of the periodic unit cube, meshed with affine hexahedra. This is separate from the first-order spherical-surface prototype; no surface higher derivatives are introduced.')
    text('Solve u<sub>t</sub> + 0.5u<sub>x</sub> + 0.3u<sub>y</sub> + 0.2u<sub>z</sub> = 0. Initially u=1 inside a ball of radius 0.2 centered at (0.30,0.35,0.40), and zero outside. At T=0.3 its exact center is (0.45,0.44,0.46). The exact field is u<sub>0</sub>(x-vT), with periodic wrapping. It varies in all three spatial directions.')
    image('volume.png',176)
    text('Point-cloud views show physical samples where u&gt;0.5, using identical cameras and axes. They locate the transported volume; transparency and this threshold hide small oscillations, which the slices below expose.','SmallPreview')
    image('slice.png',161)
    text('The same physical slice z=0.45 for both methods and the exact solution. One color range spans both numerical fields, including excursions outside [0,1]. These are polynomial evaluations, not cell averages.','SmallPreview')
    text('Configuration','HeadingPreview')
    text('Q2 tensor-product DG: 27 scalar unknowns per hexahedron; 6³ and 10³ cells (5,832 and 27,000 unknowns). MFEM Rusanov flux, classical RK4, fixed Δt=0.003, 100 steps. OFDG is applied after every full step without KXRCF or a positivity limiter. Both methods use the same nodal initialization. Initial and final numerical fields and native MFEM mesh/grid-function files are retained.','SmallPreview')
    title(19,'Oscillation reduction, with visible diffusion')
    image('cut.png',141)
    text('The right panel enlarges the leading discontinuity. Five interior samples per coordinate per cell are exported; line segments stop at cell boundaries. The cut uses exactly y=0.43, z=0.45 for both methods and the analytic reference. It is a discontinuously transported interface, not a nonlinear shock.','SmallPreview')
    tab=[['Pulse / grid','Sampled min','Sampled max','Relative L2']]
    for r in rows:
        if r['problem']==1 and r['dt']==.003:tab.append([r['method'].upper()+f" / {int(r['n'])}³",f"{r['min']:.4f}",f"{r['max']:.4f}",f"{100*r['relative_l2']:.2f}%"])
    table(tab,[180,108,108,108])
    text('<b>OFDG reduces both undershoot and overshoot but broadens the interface and increases L2 error.</b> Both pulse errors decrease on the finer grid. The method is not positivity preserving here. Error norms and extrema use order-16 tensor quadrature; integrating a discontinuous reference is approximate, and sampled extrema are not certified global extrema.','SmallPreview')
    text('Smooth and constant controls','HeadingPreview')
    text('Repeat with u<sub>0</sub>=1+0.2 sin(2πx) sin(2πy) sin(2πz). The exact solution is again translation. Errors below are absolute L2 norms, so the constant background does not dilute their interpretation.','SmallPreview')
    tab=[['Method','6³ error','10³ error','Two-grid rate']]
    for method in ['dg','ofdg']:
        a,b=[r for r in rows if r['problem']==0 and r['method']==method]
        import math
        tab.append([method.upper(),f"{a['l2']:.3g}",f"{b['l2']:.3g}",f"{math.log(a['l2']/b['l2'])/math.log(10/6):.2f}"])
    table(tab,[126,126,126,126])
    base=next(r for r in rows if r['problem']==1 and r['n']==10 and r['method']=='ofdg' and r['dt']==.003)
    half=next(r for r in rows if r['dt']==.0015)
    const=next(r for r in rows if r['problem']==2)
    text(f"Constant OFDG control: final L2 error {const['l2']:.2g}. Largest absolute periodic mass drift: {max(abs(r['mass_drift']) for r in rows):.2g}, relative to each discrete initial mass. Halving the fine pulse time step changes OFDG L2 error by {100*abs(half['l2']/base['l2']-1):.3g}%. These are small checks, not an asymptotic convergence or performance study.",'SmallPreview')
    text(f"<b>Scope and resources:</b> affine hexahedral scalar transport only; this showcase does not validate curved 3D, tetrahedra or 3D Euler. {len(run['checks'])} runs, {run['compute_seconds']:.1f} s monitored compute; peak sampled process-tree memory {run['peak_mib']:.1f} MiB. One MPI rank, one library thread; sequential runs. All checks passed. Reproduction: <b>docs/volume-showcase.md</b>; raw data and source fingerprints: <b>measurements/volume-showcase/</b>.",'SmallPreview')
