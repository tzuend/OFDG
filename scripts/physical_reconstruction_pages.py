"""Append a matched-input physical reconstruction comparison; keep earlier pages."""
import json
from reportlab.platypus import Paragraph,Spacer,Image,Table,TableStyle,PageBreak
from reportlab.lib import colors

def append_physical_pages(story,directory,styles):
 s=json.loads((directory/'summary.json').read_text());run=s['run'];c=s['controls']
 def text(t,style='SmallPreview'):story.append(Paragraph(t,styles[style]))
 def title(n,t):story.append(PageBreak());text(f'NUMERICAL PREVIEW / {n:02d} / PHYSICAL RECONSTRUCTION');text(t,'TitlePreview')
 def heading(t):text(t,'HeadingPreview')
 def figure(name):
  im=Image(str(directory/'figures'/f'{name}.png'));im.drawHeight*=504/im.drawWidth;im.drawWidth=504;story.extend([im,Spacer(1,5)])
 def table(rows,widths):
  tb=Table([[Paragraph(str(v),styles['TablePreview']) for v in row] for row in rows],colWidths=widths,hAlign='LEFT')
  tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#007f86'))]));story.extend([tb,Spacer(1,6)])
 def rep(mapping,family,order,method,p=3):return next(r for r in s['representatives'] if r['map']==mapping and r['family']==family and r['order']==order and r['method']==method and r['p']==p)
 title(12,'One physical reconstruction: does it help?')
 text('Compare the current recursive projected derivatives with exact derivatives of one physical-polynomial reconstruction. Both feasible methods receive the same DG coefficient vector. Existing OFDG code and the previous numerical studies are retained.','BodyPreview')
 heading('Construction and fair comparison')
 text('<b>Current:</b> repeatedly apply the production matrices D<sub>K,j</sub>=M<sub>K</sub><super>-1</super>G<sub>K,j</sub> in the mapped DG space. <b>Alternative:</b> reconstruct w=P<sub>q</sub>u<sub>h</sub> by physical L2 projection into total-degree physical polynomials, q=k,k+1,k+2; evaluate partial<super>alpha</super>w analytically. No intermediate projection and no permutation dependence remain. The solution and its zero-order jump need not be replaced.')
 text('The physical basis uses products of Legendre polynomials in centred, scaled physical x and y, retaining indices with total degree at most q. Weighted QR solves the projection. Quadrature uses the curved element map, but the polynomial basis is defined in physical coordinates. Derivatives of this reconstructed polynomial are exact up to roundoff, not exact derivatives of the original u<sub>h</sub>.')
 text('<b>Coverage:</b> both earlier geometry maps; P/Q2-P/Q4; all canonical derivatives through order min(k,3); three smooth fields and constants on the coupled map, plus the exactly representable quadratic on the separable map. Strong-map refinement uses n=4,8,16 for quads and n=2,4,8 for triangles; mild maps, affine controls, basis changes and quadrature sensitivity are retained in the manifest.')
 figure('ratios')
 text('Strong curvature, n=8, k=3. Values are ratios of the largest relative L2 error over derivatives of each order: physical reconstruction / current recursion. Blue is better, red worse, and 1 means equal. The worst direction can differ between methods; individual paired-direction errors are retained in comparisons.csv. Ratios avoid hiding near-coincident curves.')
 text('This is a comparison of different auxiliary spaces, not merely a basis change: physical P<sub>q</sub> has (q+1)(q+2)/2 coefficients. For Q3 input there are 16 mapped coefficients versus 10,15,21 physical reconstruction coefficients at q=3,4,5. The physical degree-k reconstruction can discard representable input content even on an affine quadrilateral.')
 title(13,'Absolute accuracy and the input limitation')
 heading('Worst-direction relative errors at n=8, k=3')
 rows=[['Field / element / order','Current','q=k','q=k+1','q=k+2','Analytic oracle']]
 for mapping,label in [('coupled','Sine'),('separable','Exact input')]:
  for family in ['quad','tri']:
   for order in [1,2,3]:
    values=[rep(mapping,family,order,m)['relative_l2'] for m in ['recursive','physical_k','physical_k1','physical_k2','oracle_k2']]
    rows.append([label+' / '+('Q3' if family=='quad' else 'P3')+' / '+str(order)]+[f'{100*v:.3g}%' for v in values])
 table(rows,[154,70,70,70,70,70])
 text('Strong coupled amplitude a=1.2 for the sine, strong separable amplitude a=0.7 for the exactly represented input. Relative errors divide each derivative L2 error by its analytic derivative norm; the table takes the largest ratio at each order. Absolute L2 errors, sampled maxima and individual derivative directions are retained in the raw results and representatives.csv.')
 heading('The oracle is diagnostic, not an available method')
 text('The last column projects the analytic solution directly into physical P<sub>k+2</sub> before differentiating. It does not receive the DG approximation and is excluded from feasible-method improvement counts. For the sine it can be far more accurate, showing that the enriched space has useful approximation capacity which reconstruction of the existing u<sub>h</sub> does not recover. In the exactly represented-input study the oracle and feasible q=k+2 reconstruction agree up to numerical error.')
 heading('Why eliminating repeated projection is not enough')
 text('The alternative differentiates one reconstruction error: partial<super>alpha</super>(P<sub>q</sub>u<sub>h</sub>-u<sub>h</sub>). It removes the intermediate projection errors but adds this initial change of approximation space. Moreover, a derivative of total order l of a degree-q physical polynomial has degree at most q-l. With k=3, q=5 and l=3, the alternative third derivative is only quadratic; the current method retains a full mapped order-3 derivative state. Enrichment by two is therefore not an equal derivative-space comparison.')
 text('<b>Observed outcome:</b> q=k can improve some smooth-field errors and worsen others. Enrichment tends to recover the current result for the projected sine rather than greatly outperform it. On the exactly represented strongly distorted input, all three tested physical reconstruction degrees have larger worst-order errors than current recursion in the displayed k=3 cases. Commuting derivatives alone are not a sufficient accuracy criterion.')
 title(14,'Face behaviour, controls and decision')
 heading('What changes at faces?')
 rows=[['Strong n=8 cases','q','Volume ratio median','Face ratio median','Jump RMS ratio median']]
 for a in s['aggregate']:
  rows.append([('Smooth fields' if a['map']=='coupled' else 'Exact input'),{'physical_k':'k','physical_k1':'k+1','physical_k2':'k+2'}[a['method']],f"{a['error_ratio_median']:.3g}",f"{a['face_ratio_median']:.3g}",f"{a['jump_ratio_median']:.3g}"])
 table(rows,[129,39,112,112,112])
 text('Ratios are alternative/current for matched derivative directions, P/Q2-P/Q4, both element families; only nonzero, resolved reference errors enter the summaries. Face error uses both one-sided traces. Jump RMS integrates the difference of the two numerical derivative traces over interior faces, where the analytic jump is zero for these smooth fields. It is a sensor-relevant diagnostic, not the full OFDG rate: no wave-speed, amplitude, factorial or cell-size weights are applied. Smaller jumps alone could also reflect lost information.')
 heading('Numerical controls and the retained failed attempt')
 text(f"Independent affine quadratic reproduction gives at most {c['affine_quadratic_max_l2']:.2g} absolute L2 derivative error. Affine triangle reconstructions agree with the mapped baseline within {c['affine_identity_max_scaled_norm_change']:.2g} in scaled error norms. Recomputed baseline norms agree with the preserved earlier study within {c['old_baseline_max_scaled_norm_change']:.2g} across {c['old_baseline_count']} matches. Reconstruction cell-integral defects are below {c['mean_defect_max']:.2g}.")
 text(f"Changing input basis (Gauss-Lobatto, Gauss-Legendre, Bernstein) changes scaled error norms by at most {c['basis_max_scaled_norm_change']:.2g}; increasing reconstruction/evaluation quadrature by eight changes them by {c['quadrature_max_scaled_norm_change']:.2g}. Nested reconstruction errors do not increase beyond {c['nested_reconstruction_violation']:.2g}. These are norm comparisons, not proof of pointwise equality. Constant-field sampled high derivatives reach {c['constant_max_sampled']:.2g}, reflecting amplified roundoff.")
 text('The initial mass-matrix solve failed the affine polynomial control: triangular degree-6 reconstruction produced a spurious third derivative of about 2.3e-5 in L2. The weighted evaluation matrix is now solved by twice-reorthogonalized QR, avoiding the conditioning penalty of normal equations. All first-attempt data and its source are retained in measurements/physical-reconstruction-normal-equations/. Tables above use the QR rerun.')
 heading('Decision for this project')
 text('<b>Keep the current recursive projection as the production approach.</b> A physical reconstruction is mathematically coherent and removes ordering dependence, but the tested degrees do not establish a consistent accuracy advantage. In particular, our exactly represented-input test favours the current approach. This is evidence against replacing it on the strength of commutation alone, not a rejection of physical-polynomial DG.')
 text('The next targeted comparison would use a reconstruction degree high enough to account for derivative degree loss (for example q=k+3 when third derivatives are required), or a derivative-aware reconstruction objective. A full physical-space DG method is a different experiment because it changes the evolved solution itself. Neither approach has been compared here in a PDE evolution; no new stability or shock-robustness claim follows.')
 text(f"<b>Resources:</b> {len(run['checks'])} QR checks, {sum(x['returncode']!=0 for x in run['checks'])} failed runs; one rank/thread on 8 logical processors. {run['compute_seconds']/60:.2f} minutes charged compute including the first attempt; QR sampled peak RSS {run['peak_mib']:.1f} MiB. Both attempts enforce 5 minutes/check and 4 GB. Commands, sources, controls and raw data: <b>measurements/physical-reconstruction/</b>. Reproduction: <b>docs/physical-reconstruction.md</b>. Earlier report pages and datasets are preserved; production code is unchanged.")
