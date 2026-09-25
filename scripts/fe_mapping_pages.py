"""Append three pages of measured FE-geometry derivative evidence."""
import json
from reportlab.platypus import Paragraph,Spacer,PageBreak,Image,Table,TableStyle
from reportlab.lib import colors

def append_fe_mapping_pages(story,directory,styles):
 summary=json.loads((directory/'summary.json').read_text());run=json.loads((directory/'run.json').read_text());highrun=json.loads((directory/'high-order/run.json').read_text());c=summary['controls'];rep=summary['representative_worst_third_relative']
 def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
 def title(number,s):story.append(PageBreak());text(f'FE GEOMETRY EXPERIMENT / {number:02d}','SmallPreview');text(s,'TitlePreview')
 def heading(s):text(s,'HeadingPreview')
 def image(name,height):story.extend([Image(str(directory/'figures'/name),width=504,height=height),Spacer(1,5)])
 def table(rows,widths):
  t=Table([[Paragraph(str(x),styles['TablePreview']) for x in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')]) ]));story.extend([t,Spacer(1,6)])
 title(20,'The mesh itself supplies<br/>the mapping derivatives')
 text('<b>Question:</b> can finite element coordinate fields replace the analytic mapping formulas needed by the earlier direct-derivative experiment? <b>Result:</b> yes for the tested polynomial volume meshes through third order. The new method agrees with analytic-map differentiation within the prescribed scaled tolerance.')
 image('meshes.png',173)
 text('Left: the coupled square map has genuinely curved element edges. Right: x=s+0.2t(1-t), y=t+0.15r(1-r), z=r, with determinant one and curved faces. Geometry is represented exactly in the selected FE spaces in exact arithmetic; the plots show the defining polynomial meshes.','SmallPreview')
 heading('A shared evaluator, with independent geometry and solution spaces')
 text('Read the geometry basis and coordinate coefficients from MFEM (GetFE and GetPointMat). Treat x, y and z as components of an FE field. Evaluate their reference-coordinate Taylor expansions with the same polynomial machinery as the solution; invert the local geometry series and compose the original solution with it. No derivative is projected along this new path. Existing geometry coefficients are reused, so no mesh projection is introduced here.')
 table([['Method','Source of geometry / derivative construction'],['Recursive','Unchanged production projected physical derivative matrices; repeat in canonical x, y, z order.'],['Analytic direct','Original mapped solution differentiated using explicit fixture mapping formulas.'],['FE geometry','Same mapped solution differentiated using only stored geometry coefficients and basis.']],[94,410])
 text('The 2D suite repeats all 50 earlier configurations. Nine small 3D cases compare Q2/Q3 solutions on 2³ and 4³ cells with Q2 geometry, including affine and quadrature controls. All methods receive identical input coefficients and evaluation points. Geometry degree and solution degree are independent.','SmallPreview')
 text('<b>Scope:</b> an isolated benchmark, not an OFDG or PDE modification. Production sources, existing benchmarks and thesis content are preserved. Agreement validates the representation-based route, not improved PDE accuracy or arbitrary geometry support.','SmallPreview')

 title(21,'Three methods on the same<br/>represented mesh and solution')
 text('Entries are the <b>worst relative volume L2 error across third derivatives</b>. The maximizing derivative may differ between methods. 2D: n=8, P/Q3, strong curvature. The exactly represented quadratic removes input approximation error in exact arithmetic; the sine is projected into the same solution space.','SmallPreview')
 data=[['Input / space','Recursive','Analytic direct','FE geometry']]
 for family,label in [('quad','Q3'),('tri','P3')]:
  for field,name in [('represented_quadratic','Exact quadratic'),('sine','Projected sine')]:
   data.append([name+' / '+label,*[f'{rep[f"{family}-{field}-{method}"]:.3g}' for method in ['recursive','direct','fe_geometry']]])
 table(data,[180,108,108,108]);image('comparison-2d.png',179)
 text('Values are fractions, not percentages. Plots show Q3 only. For exactly represented inputs the direct methods remove intermediate projection error, but retain amplified floating-point and basis-conversion error. For the sine, differentiated input approximation dominates: removing operator projection does not guarantee lower error against the analytic field. Overlapping direct-method curves are expected.','SmallPreview')
 heading('Independent checks and preserved baselines')
 table([['Check','Largest observed discrepancy'],['Earlier baseline rows',f'{c["baseline_rows"]} matched; {c["baseline_max_scaled_difference"]:.1g} change in tested norms'],['FE vs analytic direct derivatives',f'{c["method_agreement"]:.2g} on element-scaled derivative magnitudes'],['Geometry derivatives / shared faces',f'{c["geometry_derivative_residual"]:.2g} / {c["face_geometry_residual"]:.2g}'],['Independent physical finite differences',f'{c["independent_probe_count"]} probes; scaled error {c["independent_max_scaled_error"]:.2g}'],['Finite-difference step sensitivity',f'{c["independent_max_step_sensitivity"]:.2g} (steps 1e-5 and 1e-6)']],[250,254])
 text(f'Positions, Jacobians, inverse residuals, constants, affine polynomials, face traces, basis changes and quadrature checks passed. Constant derivative L2 errors reach {c["constant_fe_max_l2"]:.2g}; there is no exact-zero claim. One-sided face errors and interior jumps are retained in the data, but are not complete OFDG damping rates.','SmallPreview')

 title(22,'Curved 3D works;<br/>high orders expose numerical limits')
 table([['3D Q3, n=4: third derivatives','Recursive','Analytic direct','FE geometry'],['Exactly represented polynomial',*[f'{rep[f"hex-represented_polynomial-{m}"]:.3g}' for m in ['recursive','direct','fe_geometry']]],['Projected sine',*[f'{rep[f"hex-sine-{m}"]:.3g}' for m in ['recursive','direct','fe_geometry']]]],[210,98,98,98])
 text('Worst relative volume L2 error among nonzero third derivatives. All 19 derivative multi-indices through order three are evaluated; zero-reference derivatives use absolute errors. The exact pullback is s²+t²+r²+str. The sine is sin(2x+3y+z).','SmallPreview')
 image('high-order.png',174)
 text('<b>Higher-order stress test:</b> Q5 nodal inputs on representative elements of size h=1/2, 1/8 and 1/32; three points, all pure and mixed derivatives through order ten, affine and curved controls. References use 90-digit explicit-inverse formulas/polynomials, checked at 120 digits. Curves show max |error|/max(1, |reference|) across points and directions; this is not uniformly a relative error.','SmallPreview')
 data=[['Curved probe','h','Order 5','Order 10']]
 for dim in [2,3]:
  for n in [2,32]:
   vals=[]
   for order in [5,10]:
    a=next(r for r in summary['high_order_summary'] if r['dimension']==dim and r['n']==n and r['curved']==1 and r['field']=='quintic' and r['method']=='fe_geometry' and r['order']==order);vals.append(f'{a["max_scaled_error"]:.3g}')
   data.append([f'{dim}D',f'1/{n}',*vals])
 table(data,[126,126,126,126])
 text('<b>Interpretation:</b> the construction extends to any finite order; accurate evaluation in ordinary precision does not. Differentiation amplifies nodal initialization, polynomial conversion, geometry and arithmetic errors, particularly on small elements. Fifth-order accuracy is useful in several probes but is not uniformly near roundoff; tenth order can fail badly. The Q5 input is exactly representable mathematically, not immune to floating-point errors. Analytic-map direct evaluation shows essentially the same 3D h=1/32 errors (0.00151 at order five; 3.24e8 at order ten).','SmallPreview')
 text(f'Order-five and order-ten specializations agree through order five within {c["order5_vs_order10_max_scaled_difference"]:.2g} on a max(1, magnitude) scale. The smaller specialization limits work for ordinary use. These probes identify accuracy limits; they do not establish production performance or high-order stability.','SmallPreview')
 text(f'<b>Reproduce:</b> docs/fe-mapping-derivatives.md. {len(run["checks"])} main cases plus {len(highrun["checks"])} high-order cases; {summary["row_count"]:,} main and {summary["high_order_row_count"]:,} high-order rows. Monitored compute {run["compute_seconds"]+highrun["compute_seconds"]:.1f} s; peak sampled memory {max(run["peak_mib"],highrun["peak_mib"]):.1f} MiB. One rank/thread, sequential measured cases; independent reference analysis and compilation are additional.','SmallPreview')
