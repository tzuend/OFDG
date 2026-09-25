"""Optional two-page operator-study addendum to the existing supervisor brief."""
import json
from pathlib import Path
from reportlab.platypus import Paragraph,Spacer,Image,Table,TableStyle,PageBreak
from reportlab.lib import colors

def append_derivative_pages(story,directory,styles):
    s=json.loads((directory/'summary.json').read_text());run=s['runs']
    def text(t,style='SmallPreview'):story.append(Paragraph(t,styles[style]))
    def title(number,t):
        story.append(PageBreak());text(f'NUMERICAL PREVIEW / {number:02d} / OPERATOR ADDENDUM');text(t,'TitlePreview')
    def heading(t):text(t,'HeadingPreview')
    def figure(name):
        im=Image(str(directory/'figures'/f'{name}.png'));im.drawHeight*=504/im.drawWidth;im.drawWidth=504;story.extend([im,Spacer(1,5)])
    def table(rows,widths):
        tb=Table([[Paragraph(str(v),styles['TablePreview']) for v in row] for row in rows],colWidths=widths,hAlign='LEFT')
        tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#007f86'))]));story.extend([tb,Spacer(1,6)])
    title(7,'Curved derivatives: does order matter?')
    text('Isolated manufactured-derivative experiment; no PDE evolution and no change to the production OFDG algorithm. Exact mixed partial derivatives commute. This check measures the effect of inserting a physical L2 projection after each differentiation.','BodyPreview')
    heading('Construction and comparisons')
    text('<b>Domain and geometry:</b> start from (s,t) in [0,1]<super>2</super>; set b = a s(1-s)t(1-t), x = s+b and y = t+0.7b. Use exact degree-4 polynomial geometry with a = 0, 0.4, 1.2 (affine, mild and stronger curvature). The domain boundary stays fixed; det J is at least 1-0.425a &gt; 0.')
    text('<b>Spaces and data:</b> triangle P3-P5 and quadrilateral Q3-Q5; n = 2,4,8, plus n = 16 at a = 1.2. Triangles have 2n<super>2</super> elements, quadrilaterals n<super>2</super>; these are not equal-DOF comparisons. Physically L2-project sin(2x+3y), exp(x+y/2), and x<super>2</super>y+xy<super>2</super> into each space before differentiating.')
    text('<b>Sequences:</b> words denote chronological operations: xxy means x, then x, then y. Test xy/yx; xxy/xyx/yxx; and xyy/yxy/yyx. The canonical sequence matches the production derivative graph. Average all distinct permutations as functions, then measure physical L2 errors against the analytic derivative. The average is a candidate approximation, not a reference solution.')
    figure('convergence')
    text('<b>Refinement:</b> sine field, a = 1.2, xxy derivative. The exact derivative is -12 cos(2x+3y). Solid circles show the production order; dashed crosses show the permutation mean. Near-overlap is real; the next page plots changes directly. Both approaches improve under refinement in all six configurations.', 'SmallPreview')
    text('Errors include differentiation of the initial solution projection error as well as the intermediate derivative projections. A separate direct projection of the analytic derivative supplies an approximation floor in the raw data. This is not a comparison with the exact higher derivative of the represented DG function.', 'SmallPreview')
    title(8,'Averaging: a change, not a guaranteed gain')
    figure('averaging')
    text('Each point is one curved n=8 case: P/Q3-P/Q5, both nonzero curvature amplitudes, and xxy or xyy. Horizontal position measures the L2 change to the derivative, divided by its canonical error. Vertical position measures the resulting error ratio. Points above 1 became less accurate after averaging. Axes use different scales: the triangle accuracy changes are much smaller.', 'SmallPreview')
    rows=[['Sine, a=1.2, n=8','Canonical error','Mean error','Error change']]
    for r in s['representative']:
        rows.append([('Q' if r['family']=='quad' else 'P')+str(int(r['p']))+' / xxy',f"{r['canonical_error']:.4g}",f"{r['mean_error']:.4g}",f"{100*(r['mean_to_canonical']-1):+.2f}%"])
    table(rows,[177,109,109,109])
    text(f"Across {s['fine_third_count']} curved n=8 third-derivative cases, averaging lowers the error in {s['fine_better']}. Mean/canonical error ratios range from {s['fine_ratio_min']:.3f} to {s['fine_ratio_max']:.3f} (median {s['fine_ratio_median']:.3f}). The derivative change itself ranges from {100*s['fine_change_min']:.1f}% to {100*s['fine_change_max']:.1f}% of canonical error.")
    text('The largest worsening is mild-curvature Q3, cubic field, xyy: error 0.00109 to 0.00242 (exact derivative = 2). The relative errors are only 0.055% and 0.121%; the factor alone exaggerates its absolute significance.')
    heading('Controls and interpretation')
    text(f"Affine canonical-to-mean differences are at most {s['affine_max_change']:.2g} in L2; the affine cubic derivative error is at most {s['affine_cubic_error']:.2g}. Independent first-derivative weak residuals are below {s['first_weak_max']:.2g}; constant-derivative coefficients below {s['constant_max']:.2g}. Raising assembly quadrature by eight on four strong-curvature cases changes reported error norms by at most {s['quadrature_max_relative_error_change']:.2g} relative.")
    text('<b>Decision:</b> retain the canonical implementation for now. Averaging removes the arbitrary ordering convention, but does not reliably reduce derivative error. These volume norms do not measure OFDG face-jump rates, damping, or final PDE accuracy. A targeted face-sensor and smooth-advection comparison is the next step before changing the method.')
    text(f"<b>Resources and evidence:</b> {len(run['checks'])} sequential checks, one rank/thread; {run['compute_seconds']/60:.2f} minutes total, sampled peak {run['peak_mib']:.1f} MiB. No new large PDE runs. Raw data, commands, hashes, quadrature controls and derived tables: <b>measurements/derivative-permutations/</b>. Reproduction: <b>docs/derivative-permutations.md</b>. This addendum is separate from the thesis.")
