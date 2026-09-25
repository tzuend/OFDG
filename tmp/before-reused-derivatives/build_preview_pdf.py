#!/usr/bin/env python3
"""Build the discussion brief, including an available isolated operator addendum."""
import argparse
import json
from pathlib import Path
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(__file__).resolve().parents[1]
LABELS = {'dg': 'DG', 'ofdg': 'OFDG', 'ofdg-kxrcf': 'OFDG-KXRCF', 'oedg': 'OEDG'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--derivative-input', type=Path, default=ROOT / 'measurements/derivative-permutations')
    args = parser.parse_args()
    directory = args.input
    state = json.loads((directory / 'results.json').read_text())
    analysis = json.loads((directory / 'analysis.json').read_text())
    destination = ROOT / 'output/pdf/supervisor-preview.pdf'
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Embed the text fonts so PDF readers do not substitute a condensed font.
    regular = Path('/System/Library/Fonts/Supplemental/Arial.ttf')
    bold = regular.with_name('Arial Bold.ttf')
    if not regular.exists():
        regular = Path('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
        bold = regular.with_name('DejaVuSans-Bold.ttf')
    pdfmetrics.registerFont(TTFont('Preview', str(regular)))
    pdfmetrics.registerFont(TTFont('Preview-Bold', str(bold)))
    pdfmetrics.registerFontFamily('Preview', normal='Preview', bold='Preview-Bold')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='BodyPreview', fontName='Preview', fontSize=9.8,
                             leading=13.1, textColor=colors.HexColor('#243447'), spaceAfter=7))
    styles.add(ParagraphStyle(name='SmallPreview', parent=styles['BodyPreview'],
                             fontSize=8.7, leading=11.3, spaceAfter=6))
    styles.add(ParagraphStyle(name='TablePreview', parent=styles['SmallPreview'],
                             fontSize=8.5, leading=10.4, spaceAfter=0))
    styles.add(ParagraphStyle(name='TitlePreview', fontName='Preview-Bold', fontSize=22,
                             leading=25, textColor=colors.HexColor('#063f49'), spaceAfter=11))
    styles.add(ParagraphStyle(name='HeadingPreview', fontName='Preview-Bold', fontSize=12,
                             leading=15, textColor=colors.HexColor('#007f86'), spaceBefore=8, spaceAfter=6))
    story = []

    def text(content, style='BodyPreview'):
        story.append(Paragraph(content, styles[style]))

    def title(number, content):
        text(f'NUMERICAL PREVIEW / {number:02d}', 'SmallPreview')
        text(content, 'TitlePreview')

    def heading(content):
        text(content, 'HeadingPreview')

    def figure(name, width=504):
        image = Image(str(directory / 'figures' / f'{name}.png'))
        image.drawHeight *= width / image.drawWidth
        image.drawWidth = width
        story.extend([image, Spacer(1, 4)])

    def table(rows, widths):
        formatted = [[Paragraph(str(value), styles['TablePreview']) for value in row] for row in rows]
        item = Table(formatted, colWidths=widths, repeatRows=1, hAlign='LEFT')
        item.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e0eff0')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('LINEBELOW', (0, 0), (-1, 0), .7, colors.HexColor('#007f86')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f7f8')]),
        ]))
        story.extend([item, Spacer(1, 6)])

    def page():
        story.append(PageBreak())

    title(1, 'OFDG: a laptop-sized<br/>numerical preview')
    text('Initial evidence for choosing the next larger experiments', 'HeadingPreview')
    text('This brief compares accuracy, shock resolution, and positivity on small paper-inspired tests. It supports a supervisor discussion; it does not establish a universal method ranking.')
    heading('What the completed runs suggest')
    text('<b>Smooth accuracy:</b> on the finest P2 grids, KXRCF gating keeps errors close to DG; ungated OFDG has larger errors.<br/><b>Scalar shocks:</b> filtering suppresses the large sampled DG excursions. The two OFDG variants are often nearly coincident.<br/><b>Next decision:</b> use Shu-Osher wave preservation and curved geometry to distinguish the methods more carefully before scaling up.')
    table([
        ['Experiment', 'Configuration', 'Question'],
        ['Smooth advection', 'P1-P3; N=32/64/128/256; T=1.1', 'Accuracy across order'],
        ['Piecewise advection', 'P2; N=128/256; T=1.1', 'Oscillations and smearing'],
        ['Burgers', 'P2; smooth N=32/64/128, T=0.6;<br/>shock N=128/256, T=2.2', 'Nonlinear shock formation'],
        ['Lax', 'P1/192, P2/128, P3/96; T=1.3', 'Equal unknown count'],
        ['Shu-Osher', 'P2; N=200/400; T=1.8', 'Shock-wave interaction'],
        ['Straight / curved transport', 'P2; 8²/16²/32² cells; T=1', 'Geometry and accuracy'],
        ['2D Riemann', 'P2; 32²/64² cells; T=0.25', 'Wave structure and positivity'],
    ], [109, 226, 169])
    heading('Methods and notation')
    text('<b>DG</b> is discontinuous Galerkin without the oscillation filter. <b>OFDG</b> is its oscillation-free variant; <b>OEDG</b> is the oscillation-eliminating comparator. <b>KXRCF gating</b> applies OFDG only in cells flagged by a troubled-cell indicator. Pk denotes polynomial degree k; N is the 1D cell count. Squared counts specify 2D grids.', 'SmallPreview')
    text('All runs use RK4. CFL is 0.15 for scalar problems and 0.30 for Euler. OFDG filters after each full step; OEDG filters at RK stages. Every Euler method, including DG, uses the same positivity limiter. OEDG appears on selected finer affine shock cases only.', 'SmallPreview')
    heading('Compute used')
    text(f"<b>{analysis['run_count']} solver runs and {len(state['references'])} reference solves</b> completed in {analysis['elapsed_minutes']:.2f} minutes, with sampled peak process-tree memory {analysis['peak_mib']:.0f} MiB. Runs were sequential and used at most two MPI ranks. No cases failed or were dropped.")
    text('Ceilings: 4 GB, five minutes per check, 60 minutes total; no new checks after 45 minutes. Memory peaks are sampled. Build and document preparation are excluded from solver time. Single-run timings do not establish performance differences.', 'SmallPreview')

    page(); title(2, 'Smooth accuracy and curved geometry')
    text('Refinement reduces every smooth error. Higher slopes in the filtered cases on these grids do not establish a higher asymptotic order.', 'SmallPreview')
    figure('accuracy'); figure('nonlinear_curved_accuracy')
    rows = [['Method', 'P1 slope', 'P2 slope', 'P3 slope', 'Curved slope', 'P2 error / DG', 'Curved error / DG']]
    def last(case, method, degree):
        return [r for r in analysis['convergence'] if r['case']==case and r['method']==method and r['p']==degree][-1]
    for method in ['dg', 'ofdg', 'ofdg-kxrcf']:
        slopes = [last(case, method, degree)['observed_order'] for case,degree in
                  [('smooth_advection',1),('smooth_advection',2),('smooth_advection',3),('curved_advection_2d',2)]]
        ratios = [last(case,method,2)['l2']/last(case,'dg',2)['l2'] for case in ['smooth_advection','curved_advection_2d']]
        rows.append([LABELS[method], *[f'{value:.2f}' for value in slopes], *[f'{value:.3f}' for value in ratios]])
    table(rows, [108, 58, 58, 58, 70, 76, 76])
    text('Slopes use the finest two grids. Nominal spatial orders are 2, 3, and 4 for P1, P2, and P3. Error ratios use N=256 (1D P2) and 32² (curved P2). Nested markers indicate nearly coincident curves: the table quantifies the hidden DG / gated-OFDG difference.', 'SmallPreview')
    text('The straight and curved 2D panels share the unit-square transport problem: velocity (0.7, 0.3), sin²(pi(x+y)) data, and a static cubic deformation for the curved case. This is a 2021-paper-inspired adaptation.', 'SmallPreview')
    changes = ', '.join(f"{LABELS[m]} {100*v:.4f}%" for m,v in analysis['temporal_sensitivity'].items())
    text('Finest-P3 CFL halving changed L2 error by '+changes+'. No change exceeded 5%, so no extra series was needed. Maximum reported periodic scalar conservation drift: '+f"{analysis['periodic_conservation_max']:.1g}." , 'SmallPreview')

    page(); title(3, 'Scalar shocks: separate near-coincident curves')
    figure('scalar_shocks')
    text('Overviews locate the shaded close-ups. The bottom row subtracts OFDG, so small differences between the filtered methods remain visible. The black reference is exact translated transport; Burgers has no reference curve.', 'SmallPreview')
    rows = [['Sampled excursion', 'Transport<br/>under', 'Transport<br/>over', 'Burgers<br/>under', 'Burgers<br/>over']]
    def excursion(value):
        return '&lt;10<super>-6</super>' if value < 1e-6 else f'{value:.3f}'
    for method in LABELS:
        transport = next(r for r in analysis['extrema'] if r['case']=='piecewise_advection' and r['method']==method)
        burgers = next(r for r in analysis['extrema'] if r['case']=='shock_burgers' and r['method']==method)
        rows.append([LABELS[method], excursion(transport['undershoot']), excursion(transport['overshoot']), excursion(burgers['undershoot']), excursion(burgers['overshoot'])])
    table(rows, [148, 89, 89, 89, 89])
    text('Excursions are beyond initial-data bounds: transport [-1, 0.951057], Burgers [-0.5, 1.5]. Values below 10<super>-6</super> are grouped for display; full precision is retained in analysis.json. The plots use 21 samples per cell, with curves broken at DG interfaces; sampled extrema are not certified continuous extrema.', 'SmallPreview')
    diffs = analysis['scalar_differences']
    text('In the displayed windows, max |gated OFDG - OFDG| is '+f"{diffs['piecewise_advection']['ofdg-kxrcf']:.2g} for transport and {diffs['shock_burgers']['ofdg-kxrcf']:.2g} for Burgers. Near coincidence is a result here, not evidence that the curves were omitted.", 'SmallPreview')

    page(); title(4, 'Euler: shock fronts and wave detail')
    figure('euler_profiles')
    text('Blue bands mark shock details; the orange band marks the Shu-Osher wave detail. Shu-Osher uses P2 / 400 cells; its overview is cropped to [-1,3] from the [-5,5] domain. The equal-DOF Lax panel uses gated OFDG with P1/192, P2/128, P3/96: 384 scalar unknowns each. Its markers are cell means, not polynomial samples.', 'SmallPreview')
    ranges = analysis['wave_ranges']
    text('Wave-window peak-to-trough density range: '+', '.join(f'{LABELS[m]} {ranges[m]:.3f}' for m in LABELS)+f", WENO {ranges['reference']:.3f}. This local range is a first look at amplitude; it does not measure phase accuracy.", 'SmallPreview')
    rows = [['WENO reference pair', 'Mean absolute difference', 'Maximum difference']]
    for case,name in [('lax','Lax'),('shu_osher','Shu-Osher')]:
        values = analysis['reference_sensitivity'][case]
        rows.append([name+' 2048 vs 4096', f"{values['mean_absolute_density_difference']:.3g}", f"{values['max_density_difference']:.3g}"])
    table(rows, [204, 162, 138])
    text('These are density differences after averaging the finer WENO reference onto the coarser grid. Large localized differences near steep fronts limit pointwise judgments even when mean differences are small. WENO5 / SSPRK3 is independent of DG, but is not exact.', 'SmallPreview')
    euler = [r for r in state['runs'] if r['case']=='shu_osher' and r['n']==400 and r['status']=='ok']
    counts = ', '.join(LABELS[r['method']]+' '+str(r['metrics']['limited_elements']) for r in euler)
    text('Shu-Osher positivity cell-limit events: '+counts+'. These are cumulative interventions, not distinct cells. Fixed far-field boundaries permit changes in domain totals; those changes are not conservation-error estimates.', 'SmallPreview')

    page(); title(5, '2D Riemann: compare the same physical cut')
    text('P2 Euler, T=0.25; four-quadrant data with outflow boundaries. Maps show 64² physical cell means on one shared scale. Its range spans all displayed means, so no values are clipped.', 'SmallPreview')
    figure('riemann'); figure('riemann_slice')
    text('Both grids are evaluated at exactly y=0.740 by reconstructing each P2 polynomial from its tensor-node samples. Curves remain broken at cell interfaces. The two resolutions are separated to keep method comparisons readable.', 'SmallPreview')
    differences = analysis['cut_differences']
    text('The two OFDG traces nearly coincide. Max |gated OFDG - OFDG| along this cut: '+f"{differences['32']['ofdg-kxrcf']:.3g} on 32² and {differences['64']['ofdg-kxrcf']:.3g} on 64². The close-ups use polynomial traces; the maps above use cell means.", 'SmallPreview')
    rows = [['Method / cells', 'Min. density', 'Min. pressure', 'Limit events', 'Rejected steps']]
    for run in state['runs']:
        if run['case']!='riemann_2d':continue
        metrics = run['metrics']
        rows.append([LABELS[run['method']]+f" / {run['n']}²", *[f"{float(metrics[k]):.3g}" for k in ['min_density','min_pressure','limited_elements','rejected_steps']]])
    table(rows, [164, 85, 85, 85, 85])
    text('Minima are the solver-reported final values. Limit events count cumulative cell interventions. All three methods stayed admissible with the common positivity treatment; these coarse grids do not establish a resolution-independent ranking.', 'SmallPreview')

    page(); title(6, 'Decisions for the next compute allocation')
    table([
        ['Priority', 'Evidence from this preview', 'Next experiment / decision'],
        ['1', 'Shu-Osher filters look close; amplitude alone cannot establish phase or wave accuracy.', 'Refine the wave region; measure amplitude and phase against a checked reference. Repeat timings only for this selected case.'],
        ['2', 'Gating nearly matches DG smooth errors; curved errors are larger on the same grid counts.', 'Refine matched straight / curved problems, then a small curved Euler case. Separate geometry error from filtering.'],
        ['3', '2D OFDG variants differ little on the sampled cut, while DG uses more positivity interventions.', 'Refine the Riemann problem with boundary conditions and limiter settings matched across future comparisons.'],
        ['Later', 'No result here requires a giant test to interpret.', 'Defer blast waves, large 3D studies, and scaling studies until the scientific question and hardware are chosen.'],
    ], [45, 224, 235])
    heading('What remains uncertain')
    text('This is a small configured comparison, not a full paper reproduction. Fluxes, initialization, filtering cadence, and positivity treatment affect the outcome. Only smooth P3 advection has a time-step sensitivity check. Shock widths and the wave-detail range need broader refinement. Native NURBS is outside the preview; curved results use supported polynomial geometry.', 'SmallPreview')
    heading('Open the evidence')
    text('Run: <b>preview</b>, 103 successful solver runs and four references.<br/>Dataset: <b>measurements/study/preview/</b><br/>Instructions: <b>docs/supervisor-preview.md</b><br/>Configuration: <b>experiments/preview_manifest.json</b>', 'SmallPreview')
    text('Paths are relative to codex_dev. The dataset retains exact commands, logs, physical samples, reference CSVs, convergence tables, and source/executable fingerprints. Existing thesis assets are unchanged. Project revision '+state['provenance']['project'][:12]+'; MFEM '+state['provenance']['mfem'][:12]+'.', 'SmallPreview')
    heading('Research basis')
    references = [
        ('[1] J. Lu, Y. Liu, C.-W. Shu (2021). An Oscillation-Free Discontinuous Galerkin Method for Scalar Hyperbolic Conservation Laws. Scalar advection and Burgers examples.', '10.1137/20M1354192'),
        ('[2] Y. Liu, J. Lu, C.-W. Shu (2022). An Essentially Oscillation-Free Discontinuous Galerkin Method for Hyperbolic Systems. Lax, Shu-Osher, and equal-DOF comparisons.', '10.1137/21M140835X'),
        ('[3] M. Peng, Z. Sun, K. Wu (online 2024; issue 2025). OEDG: Oscillation-Eliminating Discontinuous Galerkin Method for Hyperbolic Conservation Laws. First four-quadrant Riemann setup.', '10.1090/mcom/3998'),
    ]
    for description, doi in references:
        text(description+' <link href="https://doi.org/'+doi+'" color="#007f86">DOI: '+doi+'</link>', 'SmallPreview')
        story.append(Spacer(1, 3))

    if (args.derivative_input / 'summary.json').exists():
        from derivative_preview_pages import append_derivative_pages
        append_derivative_pages(story, args.derivative_input, styles)

    accuracy_directory = ROOT / 'measurements/derivative-accuracy'
    if (accuracy_directory / 'summary.json').exists():
        from derivative_accuracy_pages import append_accuracy_pages
        append_accuracy_pages(story, accuracy_directory, styles)

    direct_directory = ROOT / 'measurements/direct-derivatives'
    if (direct_directory / 'summary.json').exists():
        from direct_derivative_pages import append_direct_pages
        append_direct_pages(story, direct_directory, styles)

    sphere_directory = ROOT / 'measurements/sphere-ofdg'
    if (sphere_directory / 'summary.json').exists():
        from sphere_ofdg_pages import append_surface_ofdg_pages
        append_surface_ofdg_pages(story, sphere_directory, styles)

    volume_directory = ROOT / 'measurements/volume-showcase'
    if (volume_directory / 'summary.json').exists():
        from volume_showcase_pages import append_volume_pages
        append_volume_pages(story, volume_directory, styles)

    fe_mapping_directory = ROOT / 'measurements/fe-mapping-derivatives'
    if (fe_mapping_directory / 'summary.json').exists():
        from fe_mapping_pages import append_fe_mapping_pages
        append_fe_mapping_pages(story, fe_mapping_directory, styles)

    def footer(canvas, doc):
        canvas.setStrokeColor(colors.HexColor('#d6e1e3'))
        canvas.line(45, 38, 550, 38)
        canvas.setFont('Preview', 8)
        canvas.setFillColor(colors.HexColor('#5e707d'))
        canvas.drawString(45, 25, 'OFDG / supervisor discussion / preliminary numerical evidence')
        canvas.drawRightString(550, 25, str(doc.page))

    document = SimpleDocTemplate(str(destination), pagesize=(595.28,841.89),
                                rightMargin=45,leftMargin=45,topMargin=34,bottomMargin=48)
    document.build(story,onFirstPage=footer,onLaterPages=footer)
    print(destination)


if __name__ == '__main__':
    main()
