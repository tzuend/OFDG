"""Report the isolated, explicitly defined first-order surface OFDG adaptation."""
import json
from reportlab.platypus import Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib import colors


def append_surface_ofdg_pages(story,directory,styles):
    summary=json.loads((directory/'summary.json').read_text());run=json.loads((directory/'run.json').read_text());controls=summary['controls']
    def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
    def heading(s):text(s,'HeadingPreview')
    def title(n,s):story.append(PageBreak());text(f'SURFACE DG / OFDG PROTOTYPE / {n:02d}','SmallPreview');text(s,'TitlePreview')
    def figure(name,height):story.extend([Image(str(directory/'figures'/name),width=504,height=height),Spacer(1,4)])
    def table(rows,widths):
        t=Table([[Paragraph(str(v),styles['TablePreview']) for v in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')])]))
        story.extend([t,Spacer(1,5)])
    title(15,'OFDG on a spherical surface')
    text('<b>This experiment evolves a surface-DG solution and applies an OFDG adaptation after each time step.</b> It replaces the FEM-only sphere example in this brief. Both the existing volume implementation and the separate surface-FEM dataset are retained unchanged.')
    heading('Rigid rotation with an exact solution')
    text('On the unit sphere Γ, solve u<sub>t</sub> + div<sub>Γ</sub>(v u)=0 with v=(-y,x,0). This velocity is tangential and surface-divergence-free. The exact solution is the initial field evaluated at (x cos t + y sin t, -x sin t + y cos t, z). Run to T=π/2, a quarter revolution.')
    table([['Test','Initial field','Purpose'],['Smooth rotation','u<sub>0</sub>=1+0.5x','Check refinement and the accuracy cost of filtering.'],['Discontinuous cap','u<sub>0</sub>=1 if x&gt;0.5, otherwise 0','Check oscillation suppression and edge smearing.'],['Constant control','u<sub>0</sub>=1','Check the surface operator and filter preserve constants.']],[116,178,210])
    heading('Surface fluxes and exact curved geometry')
    text('Start with octahedral triangles or cubed-sphere quads, then use the pointwise radial map X=F<sub>h</sub>/|F<sub>h</sub>|. All quadrature points lie on the exact sphere. Its 3×2 Jacobian gives tangential gradients through J(J<super>T</super>J)<super>-1</super>. Across each edge, one conservative upwind flux uses v dotted with the tangential outward co-normal. A single flux is added with opposite signs to the neighboring cells.')
    text('This radial lift is a custom geometric evaluator in the isolated prototype, not a new general MFEM surface interface. It avoids confusing geometry approximation with the filtering comparison. Solution spaces are P1 triangles and Q1 quads; each family has 96, 384 and 1536 scalar unknowns on the three tested levels.')
    heading('The first-order surface OFDG definition used here')
    figure('formula.png',101)
    text('Π<sub>K</sub> is the physical surface L2 projection into the original local DG space; g<sub>K,a</sub> are its three projected tangential-gradient components in fixed Cartesian axes. Bars denote surface means. β<sub>F</sub> is the sampled maximum absolute co-normal speed. Face jumps are physical edge averages. A is estimated from volume and edge samples; a relative near-constant guard disables damping below 1e-12.')
    text('The coefficients 1/2 and 3/2 and exponential damping of the nonconstant part are the existing order-one OFDG construction. Surface measures, co-normal speeds and tangential derivatives are the necessary adaptations. The componentwise absolute sum follows the Cartesian convention; no rotation-invariance claim is made. Damping preserves the physical cell mean. No positivity limiter or KXRCF gate is added.','SmallPreview')
    heading('Time integration and support boundary')
    text('Both methods use the same physical L2 initialization, upwind DG operator and RK4 steps: Δt ≤ 0.12 min(h<sub>K</sub>)/3, shortened to land on T. OFDG is applied once after each complete RK4 step. This is an isolated P1/Q1 prototype, not general surface support in production OFDG. Higher-order surface derivatives, arbitrary surfaces and KXRCF remain unimplemented.','SmallPreview')

    title(16,'Rotating cap: less oscillation,<br/>more smearing')
    figure('cap-sphere.png',174)
    text('512 triangles, 1536 unknowns. The left panel is the exact initial cap; both methods start from its same unfiltered physical L2 projection, which already has overshoots. The numerical panels show T=π/2 from the same camera with one shared color scale. Rendering averages values on small display triangles; the tests use quadrature, not image colors.','SmallPreview')
    figure('cap-cut.png',162)
    text('Equatorial cut at final time: northern one-sided traces, plotted separately within each element to preserve DG jumps. The right panel enlarges the edge at 30 degrees. The exact profile is one between 30 and 150 degrees. The cut does not contain every global extremum.','SmallPreview')
    data=[['Finest mesh / method','Sampled minimum','Sampled maximum','Relative L2 error']]
    for r in summary['finest_cap']:
        data.append([('Triangles' if r['family']=='tri' else 'Quads')+' / '+r['method'].upper(),f"{r['minimum']:.4f}",f"{r['maximum']:.4f}",f"{100*r['relative_l2']:.2f}%"])
    table(data,[186,106,106,106])
    heading('What the comparison establishes')
    text('The surface OFDG adaptation substantially reduces undershoot and overshoot on both element families, but does not eliminate them. It spreads the discontinuity and increases the L2 error. The result supports a working oscillation-damping mechanism on this surface; it is neither a positivity guarantee nor evidence that filtering improves every error measure.')
    text('Errors use the analytic rotated cap and positive physical quadrature of order 60. Raising assembly/input quadrature from 20 to 30 on the finest triangular cap changes reported L2 errors by about 0.053% for DG and 0.042% for OFDG. Integration of a discontinuity is numerical; these values are diagnostic, not certified continuous norms.','SmallPreview')

    title(17,'Smooth accuracy and validation')
    figure('smooth.png',168)
    data=[['Finest mesh / method','Relative L2 error','Observed L2 rate']]
    for r in summary['finest_smooth']:data.append([('Triangles' if r['family']=='tri' else 'Quads')+' / '+r['method'].upper(),f"{100*r['relative_l2']:.4g}%",f"{r['l2_rate']:.2f}"])
    table(data,[230,137,137])
    text('Rates compare the final two refinement levels as the mesh spacing halves. Both methods converge; ungated OFDG has a larger smooth error. Rates above two over this short sequence are not a superconvergence claim. Relative error divides by the full solution norm, including its constant offset.','SmallPreview')
    heading('Conservation, constants and sensitivity')
    text(f"Across the completed suite, total surface mass changes by at most {controls['mass_drift']:.2g}, measured against each run's discrete initial mass. The largest cell-mean change caused by one filter application is {controls['filter_mean_defect']:.2g}. Constant-field final L2 error is at most {controls['constant_error_l2']:.2g}. The constant-state semidiscrete residual is at most {controls['constant_operator_residual']:.2g}; opposite-side integrated co-normal flux factors agree to {controls['flux_mismatch']:.2g}.")
    text(f"Halving the finest triangular smooth time step changes DG's L2 error by {100*controls['time_sensitivity']['dg']:.2g}% and OFDG's by {100*controls['time_sensitivity']['ofdg']:.3g}%. Raising smooth assembly/evaluation quadrature by six changes those errors by {100*controls['smooth_quadrature_sensitivity']['dg']:.3g}% and {100*controls['smooth_quadrature_sensitivity']['ofdg']:.3g}%. These changes are small relative to the DG/OFDG difference; filter cadence still affects the time-discrete result.")
    heading('Retained unsuccessful measurement attempt')
    text('The initial triangle integration rules included negative weights. They produced an invalid negative squared-error integral in the constant control and unreliable cap error norms. Those measurements are not used here. The source and raw output remain archived in measurements/sphere-ofdg-initial-quadrature/. The rerun uses positive tensor/Duffy rules; its checks also reject negative or nonfinite squared-error integrals. No OFDG damping formula was changed to correct that measurement issue.','SmallPreview')
    heading('Decision and remaining work')
    text('<b>Surface OFDG is now demonstrated in an isolated first-order transport prototype.</b> Keep production volume OFDG unchanged. General integration would require a reviewed surface geometry/flux interface and a deliberate definition of higher-order tangential derivatives. Higher-order surface OFDG and general surface support are out of scope. Retain this first-order example only as an exploratory demonstration.')
    text(f"<b>Resources:</b> {len(run['checks'])} final checks, all passed; {run['compute_seconds']:.2f} s total monitored numerical time including the first attempt; peak sampled RSS {run['peak_mib']:.1f} MiB. One rank/thread, sequential; five-minute/check and 4 GB ceilings respected. Details: <b>docs/sphere-ofdg.md</b>. Commands, source fingerprints, initial/final samples and raw results: <b>measurements/sphere-ofdg/</b>. The earlier surface-FEM study is archived separately.",'SmallPreview')
