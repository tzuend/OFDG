"""Append the isolated exact mapped-derivative comparison to the supervisor brief."""
import csv
import json
from reportlab.platypus import Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib import colors


def append_direct_pages(story, directory, styles):
    summary=json.loads((directory/'summary.json').read_text())
    run=json.loads((directory/'run.json').read_text())
    controls=summary['controls']
    rows=list(csv.DictReader((directory/'results.csv').open()))
    def text(s,style='BodyPreview'):story.append(Paragraph(s,styles[style]))
    def heading(s):text(s,'HeadingPreview')
    def title(n,s):
        story.append(PageBreak());text(f'ISOLATED DERIVATIVE STUDY / {n:02d}','SmallPreview');text(s,'TitlePreview')
    def table(data,widths):
        item=Table([[Paragraph(str(x),styles['TablePreview']) for x in row] for row in data],colWidths=widths,repeatRows=1,hAlign='LEFT')
        item.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e0eff0')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')])]))
        story.extend([item,Spacer(1,6)])
    def selected(field,family,method,word='xxx'):
        return next(r for r in rows if r['field']==field and r['family']==family and r['method']==method and r['n']=='8' and r['p']=='3' and r['basis']=='gll' and r['extra']=='0' and r['derivative']==word and float(r['amplitude'])==(1.2 if field=='sine' else .7))

    title(12,'Differentiate the original<br/>mapped DG function')
    text('<b>Question:</b> can we evaluate physical derivatives of the existing DG function directly, avoiding every intermediate projection? This is an isolated operator experiment. The production OFDG implementation and all earlier raw datasets remain unchanged.')
    heading('What changes, and what stays the same?')
    table([['Current recursive projection','Direct mapped differentiation'],['After each physical derivative, project back into the original mapped DG space. Repeat on the projected coefficients.','Differentiate the original mapped function through the curved coordinate map. Evaluate the requested derivative directly at quadrature points.'],['Intermediate projections can change higher derivatives and make different derivative orderings disagree.','No derivative projection or new solution space. Mixed partials commute mathematically within each smooth, invertibly mapped element.']], [252,252])
    text('Both methods receive exactly the same physical L2-initialized DG coefficients. Direct differentiation does <b>not</b> reconstruct a polynomial in physical coordinates and does not use analytic solution derivatives to calculate the numerical answer. Analytic derivatives are used only as error references.')
    heading('Mathematical construction')
    text('Write the reference solution as u-hat(xi) = sum c<sub>j</sub> phi-hat<sub>j</sub>(xi), and the physical function as u<sub>h</sub>(x) = u-hat(F<super>-1</super>(x)). With B = J<super>-1</super>, physical differentiation in direction i is L<sub>i</sub> = sum<sub>a</sub> B<sub>ai</sub> partial/partial xi<sub>a</sub>. Repeated application must differentiate B as well as the reference polynomial.')
    text('Equivalently, at a point x<sub>q</sub> form the third-order Taylor expansion of F<super>-1</super>(x<sub>q</sub> + h z). Compose the original reference polynomial with that expansion. Its coefficient of z<sub>x</sub><super>a</super>z<sub>y</sub><super>b</super>, multiplied by a! b! / h<super>a+b</super>, gives the physical derivative with a x-differentiations and b y-differentiations. Truncation to degree three retains all derivatives through degree three; this is not a finite-difference approximation.')
    heading('Implementation and support boundary')
    text('The prototype uses truncated bivariate Taylor arithmetic and four formal inverse-map corrections with the base-point Jacobian. Reference polynomials are recovered in a local monomial basis by interpolation, an algebraic change of basis checked against MFEM values and first derivatives. The known polynomial geometry formula is used directly and checked against MFEM point evaluations. No global inverse-map formula is required.')
    text('Coverage: polynomial 2D triangles and quadrilaterals, P/Q2-P/Q4, derivatives through min(k,3), and Gauss-Lobatto, Gauss-Legendre and Bernstein coefficient bases. These tests do not supply a general MFEM higher-geometry-derivative interface, NURBS/3D support, a cached production implementation, or a PDE stability result. Different derivative permutations are implicit in one Taylor coefficient; separate order-by-order implementations were not compared.','SmallPreview')

    title(13,'What improves, and what does not')
    heading('Matched-input third derivatives: n=8, solution order k=3')
    data=[['Input / space','Recursive relative error','Direct relative error','Recursive / direct absolute error']]
    for field,label in [('represented_quadratic','Exact input'),('sine','Projected sine')]:
        for family in ['quad','tri']:
            a=selected(field,family,'recursive');b=selected(field,family,'direct')
            data.append([label+' / '+('Q3' if family=='quad' else 'P3'),f"{100*float(a['relative_l2']):.4g}%",f"{100*float(b['relative_l2']):.4g}%",f"{float(a['error_l2']):.3g} / {float(b['error_l2']):.3g}"])
    table(data,[125,115,115,149])
    text('Every table entry refers to the same xxx derivative. Relative error is ||numerical derivative - analytic derivative||<sub>L2(K summed)</sub> / ||analytic derivative||<sub>L2(K summed)</sub>. Absolute errors are physical volume L2 norms. Exact-input third derivatives can be large; percentages and absolute errors therefore answer different questions.','SmallPreview')
    image=Image(str(directory/'figures/comparison.png'),width=504,height=196);story.extend([image,Spacer(1,5)])
    text('Plots use the largest relative error over xxx, xxy, xyy, yyy; the maximizing direction can differ between methods. Left: direct errors are limited by amplified floating-point/input roundoff, so their increase under refinement is not a truncation-error convergence rate. Right: ratios above one mean direct differentiation has a larger error. Ratios expose differences that would be hidden by overlapping error curves.','SmallPreview')
    heading('Why these outcomes coexist')
    text('Direct differentiation removes the recursive operator error relative to the represented u<sub>h</sub>. It does not remove the differentiated input error, partial<super>alpha</super>(u<sub>h</sub>-u). For the exactly represented input, that term vanishes in exact arithmetic and the improvement is dramatic. For the projected sine, differentiation of the existing DG approximation dominates; removing intermediate projections gives little benefit and can slightly increase error. Projection can also suppress components of the input error. No general ranking against the analytic solution follows.')
    heading('Diagnostic fields and geometry')
    text('Original coordinates (s,t) cover the unit square, with n<super>2</super> quads or 2n<super>2</super> triangles. Sine: u=sin(2x+3y), x=s+b, y=t+0.7b, b=1.2s(1-s)t(1-t). Exact input: u=s(x)<super>2</super>+s(x)t(y)+t(y)<super>2</super>, x=s+0.7s(1-s), y=t+0.49t(1-t). The latter is a non-affine map with straight quad edges. Both maps are represented exactly by degree-four geometry, up to floating-point error.','SmallPreview')

    title(14,'Face diagnostics, verification<br/>and next decision')
    heading('Sensor-relevant changes')
    data=[['Strong maps, n=8','Volume error ratio','Face error ratio','Interior jump ratio']]
    for a in summary['aggregate']:
        data.append(['Smooth projected fields' if a['map']=='coupled' else 'Exactly represented input',f"{a['error_l2_median_ratio']:.4g}",f"{a['face_rms_median_ratio']:.4g}",f"{a['jump_rms_median_ratio']:.4g}"])
    table(data,[177,109,109,109])
    text('Entries are medians of direct/recursive ratios for matched derivative directions: 138 comparisons for sine, shortwave sine and exponential inputs; 46 for the represented quadratic, across P/Q2-P/Q4. Face errors use both one-sided traces; jump RMS uses interior faces, where the analytic jumps vanish. Cases with baseline norms below 1e-9 are excluded from ratios. These are derivative diagnostics, not full OFDG damping rates.','SmallPreview')
    text('For the smooth projected inputs, the volume ratios range from 0.983 to 1.373 and the jump ratios from 0.565 to 1.996. Direct differentiation is therefore not uniformly better for those fields. For the exactly represented input, all 46 volume, face and jump comparisons improve. A smaller jump alone does not establish a better shock sensor.')
    heading('Independent checks and numerical limits')
    text(f"<b>Higher derivatives:</b> 18 first-to-third derivative probes on both maps agree with independently inverted, 65-digit physical finite differences within {controls['independent_fd_max_scaled_difference']:.2g}, scaled by max(1, reference magnitude). Changing the finite-difference step from 1e-5 to 1e-6 changes the reference by at most {controls['independent_fd_max_step_sensitivity']:.2g}. Finite differences validate the jet method; they are not used by it.")
    text(f"<b>Representation and baseline:</b> polynomial evaluation differs from MFEM by at most {controls['representation_residual']:.2g}; physical first derivatives differ by {controls['first_derivative_residual']:.2g} on a max(1, magnitude) scale. Geometry positions differ by at most {controls['geometry_residual']:.2g} after division by h. The {controls['baseline_matches']} matched recursive error norms reproduce the prior accuracy study exactly in the saved CSV values.")
    text(f"<b>Sensitivity:</b> changing input basis changes error norms by at most {controls['basis_max_scaled_norm_change']:.2g}, scaled by max(1, analytic derivative L2 norm); eight extra evaluation quadrature orders change them by {controls['quadrature_max_scaled_norm_change']:.2g}. Direct represented-input relative errors are at most {controls['represented_max_relative']:.2g} across the suite. Constants have up to {controls['constant_max_l2']:.2g} absolute derivative L2 error: differentiation amplifies initialization and basis-conversion roundoff. No exact-zero claim is made.")
    heading('Decision')
    text('<b>Retain the existing production method; keep direct differentiation as a promising prototype.</b> It resolves the projection-error mechanism decisively for exactly represented inputs, but does not generally improve derivatives of projected smooth inputs. Before adoption, compare actual OFDG damping rates and PDE behaviour, including shocks, and assess geometry-interface support, cache memory and runtime. The current experiment measures neither stability nor production speed.')
    text(f"<b>Reproduction:</b> docs/direct-derivatives.md; measurements/direct-derivatives/ contains all commands, raw rows, controls and source fingerprints. {len(run['checks'])} configurations, {summary['row_count']} rows, zero failed runs; {run['compute_seconds']:.1f} s of monitored numerical runs, one rank/thread on 8 logical processors, peak sampled RSS {run['peak_mib']:.1f} MiB. The independent probes are additional small checks. All monitored cases stayed below 5 minutes and 4 GB.",'SmallPreview')
    text('The previous physical-polynomial reconstruction study is retained separately in measurements/physical-reconstruction/ and docs/physical-reconstruction.md, with its earlier PDF archived. Its pages are omitted here, including from the appendix, to keep this discussion focused.','SmallPreview')
