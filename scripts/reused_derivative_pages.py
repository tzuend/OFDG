"""Append measured basis-reuse results without changing earlier preview pages."""
import csv,json,math
from reportlab.platypus import Paragraph,Spacer,PageBreak,Image,Table,TableStyle
from reportlab.lib import colors
METHODS=['recursive','direct','fe_geometry','reused','reused_centered']
def append_reused_pages(story,directory,styles):
 summary=json.loads((directory/'summary.json').read_text());run=json.loads((directory/'run.json').read_text());c=summary['controls'];rows=list(csv.DictReader((directory/'results.csv').open()))
 def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
 def title(n,s):story.append(PageBreak());text(f'REFERENCE-MATRIX REUSE / {n:02d}','SmallPreview');text(s,'TitlePreview')
 def table(rows,widths):
  t=Table([[Paragraph(str(x),styles['TablePreview']) for x in r] for r in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')])]));story.extend([t,Spacer(1,7)])
 def image(name,height):story.extend([Image(str(directory/'figures'/name),width=504,height=height),Spacer(1,6)])
 def number(x):return f'{float(x):.2g}'
 def high(dim,n,order,method):return next(r for r in summary['high_summary'] if (r['dimension'],r['n'],r['curved'],r['field'],r['order'],r['method'])==(str(dim),str(n),'1','quintic',str(order),method))
 title(23,'Reuse the affine machinery;<br/>differentiate by the chain rule')
 text('Two additional variants reuse the existing basis-independent reference derivative matrices for both the solution and the coordinate fields. They differentiate the inverse Jacobian algebraically, then apply the chain rule repeatedly. No analytic map or inverse-map Taylor expansion enters the new evaluator.')
 table([['Stage','Construction'],['Reference derivatives','Call the unchanged production assembler D<sub>a</sub> = V<super>−1</super>G<sub>a</sub>; repeatedly apply D<sub>a</sub> to coefficient vectors. Geometry and solution may have different bases and degrees.'],['Inverse Jacobian','Scale K = J/h and B = K<super>−1</super>. Differentiate KB = I and solve for each higher derivative of B. Preserve matrix multiplication order.'],['Physical derivatives','Apply h ∂/∂x<sub>i</sub> = Σ<sub>a</sub>B<sub>ai</sub> ∂/∂ξ<sub>a</sub>, retaining product-rule terms through the requested order. Divide by h<super>m</super> at order m.'],['Optional centering','Subtract c q before differentiation, with Vq = 1. Apply separately to each coordinate and solution component; this remains valid when q is not an all-ones vector.']],[106,398])
 text('The graph driver is adapted locally to arbitrary requested order, including structural P/Q polynomial zeros. The production matrix assembler is reused literally. This experiment makes no production-interface or damping changes.','SmallPreview')
 table([['Validation','Largest observed discrepancy'],['MFEM Jacobian / physical first derivative',f'{number(c["reused_geometry_residual"])} / {number(c["reused_first_residual"])}'],['Differentiated inverse identity KB = I',number(c['reused_inverse_residual'])],['Fifth-order polynomial basis reproduction',f'{number(c["basis_polynomial_error"])} across {c["basis_cases"]} basis/centering controls'],['Preserved older rows',f'{c["baseline_rows"]:,} main + {c["high_baseline_rows"]:,} high-order values; zero measured change']],[255,249])
 text('Basis controls include Gauss-Lobatto, Gauss-Legendre and Bernstein, plus scaled bases whose constant coefficients are not all ones. The largest constant-vector representation error is '+number(c['constant_coefficient_error'])+'.','SmallPreview')
 data=[['Third-order relative L2','Recursive','Analytic','FE series','Reuse','Centered']]
 for family,n,field,label in [('quad','8','represented_quadratic','2D Q3 exact input'),('tri','8','represented_quadratic','2D P3 exact input'),('hex','4','represented_polynomial','3D Q3 exact input'),('hex','4','sine','3D Q3 projected sine')]:
  vals=[]
  for method in METHODS:
   rs=[float(r['relative_l2']) for r in rows if r['family']==family and r['n']==n and r['p']=='3' and r['field']==field and r['basis']=='gll' and r['extra']=='0' and float(r['amplitude'])>0 and len(r['derivative'])==3 and r['method']==method];vals.append(number(max(v for v in rs if math.isfinite(v))))
  data.append([label,*vals])
 table(data,[144,72,72,72,72,72])
 text('Worst nonzero-reference third derivative across each selected group (2D includes the tested curvature amplitudes). Values are fractions. The sine is limited by its FE input approximation; direct evaluation cannot recover derivatives discarded by that approximation.','SmallPreview')

 title(24,'Fifth-order accuracy improves;<br/>tenth order remains limited')
 text('All five methods receive the same Q5 coefficients and mesh. Curved local-element probes use three points and every derivative multi-index through order ten. The ideal-field reference is computed with 90-digit explicit-inverse formulas.','SmallPreview')
 image('high-orders.png',277)
 data=[['Curved, h=1/32','Recursive','Analytic','FE series','Reuse','Centered']]
 for dim,order in [(2,5),(3,5),(2,10),(3,10)]:data.append([f'{dim}D / order {order}',*[number(high(dim,32,order,m)['max_scaled_error']) for m in METHODS]])
 table(data,[144,72,72,72,72,72])
 text('At h=1/32, order-five <b>sampled relative L2</b> falls from 3.34·10<super>−8</super> (FE series) to 3.84·10<super>−9</super> (reuse) in 2D, and from 5.65·10<super>−6</super> to 1.18·10<super>−7</super> in 3D. Reuse gives maximum h<super>5</super>-scaled absolute errors of 3.55·10<super>−14</super> and 5.63·10<super>−13</super>, respectively.','SmallPreview')
 text('Plot/table metric: max |error| / max(1, |ideal derivative|), over points and directions. This is a mixed absolute/relative scale, not uniformly a relative error. The data also contain absolute errors, relative errors, order-group discrete relative L2, and h<super>m</super>-scaled errors. No integration weights enter these local stress norms.','SmallPreview')
 text('Reusing the native derivative matrices avoids the earlier monomial-conversion path. Centering reduces cancellation in some cases, especially constants, but does not improve every derivative. Both approaches retain coefficient roundoff, differentiation conditioning and the physical h<super>−m</super> amplification. Fifth-order results must be assessed by dimension and mesh size; tenth order has no general accuracy guarantee.','SmallPreview')
 text(f'Requesting capacity five instead of ten changes the shared outputs by at most {number(c["capacity5_vs10"])} on a max(1, magnitude) scale. The order-five path therefore avoids unnecessary high-order work. These are diagnostic measurements, not an optimized production implementation.','SmallPreview')

 title(25,'Separate evaluator error<br/>from the represented input')
 text('A second, independent reference differentiates the <b>same stored FE function and geometry</b>: reconstruct their tensor Lagrange polynomials from exact binary64 coefficients/nodes, invert the map with 90-digit Newton iteration, then take high-precision physical finite differences. No benchmark chain-rule or inverse-series arithmetic is shared.')
 image('stored-errors.png',151)
 data=[['Order 5, h=1/32','Analytic','FE series','Reuse','Centered']]
 for dim in ['2','3']:
  vals=[next(r for r in summary['stored_summary'] if (r['dimension'],r['n'],r['order'],r['method'])==(dim,'32','5',m)) for m in METHODS[1:]]
  data.append([dim+'D evaluation error',*[number(r['absolute_evaluation_error']) for r in vals]])
  data.append([dim+'D representation error',*[number(r['absolute_representation_error']) for r in vals]])
 table(data,[164,85,85,85,85])
 text('Table: maximum absolute errors over the selected fifth-order directions. Representation error is shared because the input is shared. Each column can attain its maximum in a different direction; these maxima must not be added. The plot uses max |evaluation error| / max(1, |stored reference|).','SmallPreview')
 table([['Measure','What it diagnoses'],['Total error','Computed derivative minus the ideal analytic-field derivative. This is the final accuracy experienced by the user.'],['Evaluation error','Computed derivative minus the high-precision derivative of stored u<sub>h</sub> on stored F<sub>h</sub>. Tests the differentiation machinery.'],['Representation error','Stored-function derivative minus the ideal-field derivative. Here the spaces can represent the fixtures exactly, so coefficient/geometry roundoff remains.'],['Relative and element-scaled errors','Relative L2 normalizes by reference norm; use absolute error for zero truth. h<super>m</super>-scaled error exposes amplification but does not replace physical error.']],[122,382])
 text('<b>Accuracy limit:</b> selected fifth-order evaluation errors exceed 10<super>−7</super> on the normalized scale (up to 7.6·10<super>−7</super> for centered reuse). The stable reference steps identify a numerical limitation, not stencil uncertainty.','SmallPreview')
 text(f'{c["stored_probe_count"]} independent pure/mixed probes, orders 1, 3 and 5, in 2D/3D at h=1/8 and 1/32. Steps h·10<super>−4</super> and h·5·10<super>−5</super> differ by at most {number(c["stored_step_sensitivity"])} on a max(1, reference) scale. This validates selected directions, not every tenth-order output.','SmallPreview')
 text(f'<b>Reproduce:</b> docs/reused-derivatives.md. {len(run["checks"])} monitored configurations; {summary["rows"]:,} main and {summary["high_rows"]:,} high-order rows. Sequential one-rank/one-thread compute: {run["compute_seconds"]:.1f} s; peak sampled memory {run["peak_mib"]:.1f} MiB. Initialization and evaluation timings are retained for the combined diagnostic program. Reference analysis and compilation are additional.','SmallPreview')
