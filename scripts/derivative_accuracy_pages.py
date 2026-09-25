"""Three-page analytic-derivative accuracy addendum for the supervisor preview."""
import json,math
from reportlab.platypus import Paragraph,Spacer,Image,Table,TableStyle,PageBreak
from reportlab.lib import colors

def append_accuracy_pages(story,directory,styles):
 s=json.loads((directory/'summary.json').read_text());run=s['run']
 def text(t,style='SmallPreview'):story.append(Paragraph(t,styles[style]))
 def title(n,t):story.append(PageBreak());text(f'NUMERICAL PREVIEW / {n:02d} / DERIVATIVE ACCURACY');text(t,'TitlePreview')
 def heading(t):text(t,'HeadingPreview')
 def figure(name):
  im=Image(str(directory/'figures'/f'{name}.png'));im.drawHeight*=504/im.drawWidth;im.drawWidth=504;story.extend([im,Spacer(1,5)])
 def table(rows,widths):
  tb=Table([[Paragraph(str(v),styles['TablePreview']) for v in row] for row in rows],colWidths=widths,hAlign='LEFT');tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#007f86'))]));story.extend([tb,Spacer(1,6)])
 def lookup(which,family,p,order,field):return next(x for x in s[which] if x['family']==family and x['p']==p and x['order']==order and x['field']==field)
 title(9,'How accurate are the derivatives?')
 text('An isolated operator study compares the production projected-derivative matrices with analytic physical derivatives. No PDE is advanced and no production filtering code is changed. Higher derivative errors can be substantial even when first derivatives are accurate.','BodyPreview')
 heading('Two questions, tested separately')
 text('<b>A. Complete approximation:</b> physically project sin(2x+3y), sin(8x+6y), and exp(x+y/2), then differentiate recursively. Errors include the input projection and all intermediate derivative projections. Use the coupled curved-face map from page 7, with a=0, 0.4, 1.2.')
 text('<b>B. Derivative projection alone:</b> use an exactly representable mapped quadratic and an analytically invertible non-affine map (next page). This removes input approximation error and permits comparison with exact derivatives of the represented DG function itself.')
 text('<b>Coverage:</b> P/Q1-P/Q4 for A; P/Q2-P/Q4 for B. All canonical coordinate derivatives through min(p,3): x, y; xx, xy, yy; xxx, xxy, xyy, yyy. Gauss-Lobatto, Gauss-Legendre and positive Bernstein bases. Refinement uses n=4,8,16 for quads and n=2,4,8 for triangles, with selected n=16 triangle checks. Basis comparisons use n=8; the raw manifest identifies every case. Geometry degree 4 represents both maps exactly. Spaces at equal n are not equal-DOF comparisons.')
 figure('accuracy')
 text('<b>Complete-approximation errors:</b> sin(2x+3y), strong coupled curvature a=1.2, Gauss-Lobatto basis. Each curve shows the largest relative volume L2 error among derivatives of the indicated order. Different panel scales expose the refinement trend; high derivative orders remain much less accurate. These are canonical sequences, not permutation averages.')
 heading('What “orders of magnitude off” means here')
 text('Absolute error E = ||D<super>alpha</super>u<sub>h</sub> - partial<super>alpha</super>u||<sub>L2</sub>; relative error R = E / ||partial<super>alpha</super>u||<sub>L2</sub>. Thus R=10<super>-2</super> means a 1% norm error, two orders below the derivative norm. It is not a pointwise guarantee of two correct digits. Zero analytic derivatives have undefined relative error and are reported with absolute errors. No regularizing denominator hides those cases.')
 title(10,'Separate projection error from input error')
 text('For the isolation test, take x=s+a s(1-s), y=t+0.7a t(1-t), and u(x,y)=s(x)<super>2</super>+s(x)t(y)+t(y)<super>2</super>. Test a=0.3 and 0.7. On every element its pullback is quadratic, so P/Q2 and above represent it exactly. Its physical derivatives are generally nonpolynomial. The inverse map and its first three derivatives are analytic; see the reproducibility note.')
 text('The separable map stretches coordinate lines nonlinearly. Quad edges stay straight, while triangle diagonals can curve; it isolates non-affine parameterization rather than reproducing the curved quad faces of study A. Both maps are nonsingular throughout the domain.')
 figure('isolated')
 text('<b>Projection-only errors:</b> strong separable map a=0.7. Initial physical projection error is at most '+f"{s['represented_input_max']:.2g}"+' in L2 across the tested cases. All plotted errors can therefore be attributed to derivative representation, recursion and numerical roundoff, rather than an unresolved input field.')
 rows=[['n=8, strong map','First derivative','Second derivative','Third derivative']]
 for fam in ['quad','tri']:
  for p in [2,3,4]:
   cells=[('Q' if fam=='quad' else 'P')+str(p)]
   for order in [1,2,3]:
    if order>p:cells.append('not used');continue
    x=lookup('isolated',fam,p,order,'represented_quadratic');cells.append(f"{x['abs_l2']:.2e}<br/>{100*x['rel_l2']:.3g}% ({x['derivative']})")
   rows.append(cells)
 table(rows,[120,128,128,128])
 text('Each cell gives absolute L2 error, then relative error; the displayed derivative has the largest relative error at that order. The worst absolute-error direction may differ. A direct physical projection of each exact derivative supplies the best-approximation floor. Higher recursive errors can exceed this floor because earlier projection errors are differentiated again. The raw results include both the floor and the excess over it.')
 title(11,'Implications for OFDG and checks')
 heading('How large are complete-approximation errors?')
 rows=[['n=8; a=1.2','Order 1','Order 2','Order 3','Worst face error']]
 for fam in ['quad','tri']:
  for p in [2,3,4]:
   cells=[('Q' if fam=='quad' else 'P')+str(p)+' / sine'];rr=[]
   for order in [1,2,3]:
    if order>p:cells.append('not used');continue
    x=lookup('representatives',fam,p,order,'sine');rr.append(x);cells.append(f"{x['abs_l2']:.2e}<br/>{100*x['rel_l2']:.3g}%")
   cells.append(f"{100*max(x['face_rel_max'] for x in rr):.3g}%");rows.append(cells)
 table(rows,[107,99,99,99,100])
 text('Volume cells show absolute L2 error / largest relative error over the canonical derivatives at each order. The last column is the largest relative one-sided face L2 error over all tested derivative orders, normalized by the analytic trace norm. It uses every element boundary, counting interior faces twice. It is a trace diagnostic, not the OFDG jump or damping rate.')
 short=[x['rel_l2'] for x in s['representatives'] if x['field']=='shortwave'];smooth=[x['rel_l2'] for x in s['representatives'] if x['field'] in ['sine','exponential']]
 text(f"Across these n=8 P/Q2-P/Q4 strong-curvature cases, the order-wise worst relative errors range from {min(smooth):.2g} to {max(smooth):.2g} for the sine/exponential fields, and {min(short):.2g} to {max(short):.2g} for the shorter wave. These are envelopes of a configured test set, not a statistical statement about typical meshes. Sampled maximum errors and nominal-h-scaled errors are also saved.")
 heading('Basis, integration and analytic controls')
 c=s['controls']
 text(f"Changing among Gauss-Lobatto, Gauss-Legendre and Bernstein bases changes the reported L2 error norms by at most {c['basis']['max_absolute_change']:.2g} absolute ({c['basis']['max_nonconstant_scaled_change']:.2g} on nonconstant fields after scaling by the larger analytic derivative/input norm). This compares error norms, not pointwise equality of the derivative functions. Independent first-derivative weak residuals are at most {s['weak_max']:.2g}. Constant-field sampled derivative magnitudes reach {s['constant_linf_max']:.2g}: repeated differentiation amplifies tiny initialization/assembly roundoff. Relative errors for these zero derivatives are undefined.")
 text(f"Raising assembly quadrature by eight changes error norms by at most {c['assembly']['max_scaled_change']:.2g} on that scale; independently raising initialization/evaluation quadrature by eight changes them by {c['evaluation']['max_scaled_change']:.2g}. Evaluation rules otherwise use production order +16. Analytic inverse-map derivatives are checked independently by high-precision finite differences. Sampled maxima are not certified continuous suprema.")
 heading('Assessment: useful, with a resolution-dependent boundary')
 text('<b>Retain this construction for OFDG, but do not describe its higher derivatives as uniformly accurate.</b> Refinement and increased order improve the tested errors; changing coefficient basis does not cure them. First derivatives can be several orders below the true derivative norm while the highest derivative used by a low-order space remains at percent-to-tens-of-percent error. Strong distortion and shorter waves require particular care.')
 text('OFDG uses derivatives as roughness indicators, not as a replacement for the DG flux. That makes this a sensible construction to continue evaluating, but volume accuracy alone cannot establish acceptable damping or PDE accuracy. Face errors can be larger. Before claiming the extension is sufficiently accurate, compare face-jump rates and smooth-advection convergence against an independent derivative construction; shock robustness requires a separate check.')
 text(f"<b>Resources:</b> {len(run['checks'])} completed checks, {sum(c['returncode']!=0 for c in run['checks'])} failures; one rank/thread on an 8-logical-processor laptop. {run['compute_seconds']/60:.2f} minutes charged compute (including a conservative interruption allowance), peak sampled RSS {run['peak_mib']:.1f} MiB. Raw data, exact commands, controls, source hashes and tables: <b>measurements/derivative-accuracy/</b>. Instructions: <b>docs/derivative-accuracy.md</b>. No new shock/PDE runs; no claim of stability or a universal error bound.")
