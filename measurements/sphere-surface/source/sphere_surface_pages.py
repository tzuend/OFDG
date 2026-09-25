"""Append an explicitly labelled surface-FEM example to the supervisor brief."""
import json
from reportlab.platypus import Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib import colors


def append_sphere_pages(story, directory, styles):
    summary=json.loads((directory/'summary.json').read_text())
    run=json.loads((directory/'run.json').read_text())
    def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
    def heading(s):text(s,'HeadingPreview')
    def title(n,s):
        story.append(PageBreak());text(f'EMBEDDED-SURFACE EXAMPLE / {n:02d}','SmallPreview');text(s,'TitlePreview')
    def figure(name,height):story.extend([Image(str(directory/'figures'/name),width=504,height=height),Spacer(1,5)])
    def table(data,widths):
        t=Table([[Paragraph(str(v),styles['TablePreview']) for v in row] for row in data],colWidths=widths,repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')])]))
        story.extend([t,Spacer(1,5)])
    title(15,'Solving a PDE on a<br/>spherical surface')
    text('<b>Purpose:</b> demonstrate and verify a solve on a curved two-dimensional mesh embedded in three-dimensional space. This is the surface itself, not a mesh filling the sphere. The stationary screened-Poisson problem works well as an inexpensive first example with a known solution.')
    heading('Problem: -Δ<sub>Γ</sub>u + u = f on the unit sphere')
    text('Choose u(x,y,z)=xy and f=7xy on Γ={x<super>2</super>+y<super>2</super>+z<super>2</super>=1}, following the pinned MFEM Example 7. The spherical Laplacian satisfies Δ<sub>Γ</sub>(xy)=-6xy. The sphere has no boundary, so no boundary conditions are imposed; the +u term removes the constant nullspace of the pure Laplacian.')
    figure('sphere.png',207)
    text('Computed quadratic solution and signed error on 512 curved triangles (refinement level 3, geometry degree 2). Left: original element edges on the visible side. Surface patches are subdivided only for rendering; colors average the three sampled vertex values per displayed patch. The coordinates come from the actual approximate surface, without projecting the rendered sphere to unit radius. The two panels use different color scales.','SmallPreview')
    heading('How surface geometry enters the solve')
    text('The element map F has two reference coordinates and three physical coordinates. Its Jacobian J is therefore 3×2, not square. The metric G gives the tangential gradient and surface measure. We solve the weak problem below with continuous H1 finite elements:')
    figure('surface-math.png',59)
    text('MFEM supplies these surface transformations in its diffusion and mass integrators. High-order mesh nodes are snapped to the unit sphere after refinement. Between nodes, the polynomial surface Γ<sub>h</sub> approximates rather than exactly equals Γ.')
    heading('Scope and error reference')
    text('Away from the exact sphere, use the radial extensions u<super>e</super>(x,y,z)=xy/(x<super>2</super>+y<super>2</super>+z<super>2</super>) and f<super>e</super>=7u<super>e</super>. Reported errors are integrated on Γ<sub>h</sub> against u<super>e</super>; they include both discretization and geometry effects. This is a conforming surface-FEM result, not a DG/OFDG or surface-transport validation. Production OFDG remains unchanged.')

    title(16,'Accuracy and geometry checks')
    figure('convergence.png',174)
    text('p is solution degree; g is geometry degree. Triangles start from an octahedron, quads from a cube. Each refinement multiplies the element count by four. The p=g=2 series includes one additional refinement. These are not equal-DOF comparisons.','SmallPreview')
    data=[['Elements / degree','Finest elements','Unknowns','Relative L2 error','L2 rate']]
    for r in summary['finest']:
        data.append([('Triangles' if r['family']=='tri' else 'Quads')+f" / {int(r['p'])}",int(r['elements']),int(r['dofs']),f"{100*r['relative_l2']:.4g}%",f"{r['l2_rate']:.2f}"])
    table(data,[143,83,75,115,88])
    text('Rate = log<sub>2</sub>(E<sub>coarse</sub>/E<sub>fine</sub>) using absolute L2 errors on the final two levels. Observed rates approach p+1; tangential H1-seminorm rates are approximately p (retained in the data). This small manufactured example does not establish behavior for arbitrary surfaces or singular solutions.','SmallPreview')
    heading('Geometry matters: hold solution degree p=2 fixed')
    data=[['Level 3 / relative L2 error','Planar g=1','Curved g=2','Curved g=3']]
    for family in ['tri','quad']:
        vals=[next(r for r in summary['geometry_comparison'] if r['family']==family and r['g']==g) for g in [1,2,3]]
        data.append(['Triangles' if family=='tri' else 'Quads']+[f"{100*r['relative_l2']:.4g}%" for r in vals])
    table(data,[201,101,101,101])
    text('Curved geometry markedly improves this test over planar facets. Raising g from 2 to 3 reduces sphere-area error but increases the displayed solution error: geometry and solution-error contributions can interact or cancel. Better geometric accuracy alone is not a guarantee of monotonically smaller solution error.','SmallPreview')
    heading('Checks, resources and reproducibility')
    c=summary['controls']
    text(f"All {len(run['checks'])} cases passed. A constant solution u=f=1 has relative L2 error at most {c['constant_relative_l2_max']:.2g}. Raising assembly/evaluation quadrature by six changes L2 error by at most {c['quadrature_scaled_l2_change_max']:.2g}; the largest algebraic relative residual is {c['algebraic_relative_residual_max']:.2g}. The closed-surface reaction balance |∫(u<sub>h</sub>-f<super>e</super>)dS| is at most {c['mean_balance_defect_max']:.2g}.",'SmallPreview')
    text(f"Sequential one-rank, one-thread execution: {run['compute_seconds']:.2f} s monitored compute, {run['peak_mib']:.1f} MiB peak sampled RSS. Every check stayed below 5 minutes and 4 GB. Details and exact configurations: <b>docs/sphere-surface.md</b> and <b>measurements/sphere-surface/</b>. Both triangle and quad solutions are saved as MFEM mesh/grid-function files for GLVis. Adapted from the installed <b>mfem/examples/ex7.cpp</b>; no MFEM update was made.",'SmallPreview')
