# Supervisor preview: independent presentation review

Two-agent workflow: one critical judge reads and visually inspects the PDF;
one revising agent changes the reproducible analysis and document source.
The stopping rule is a judge score of at least 8.5/10, with at most three reviews.
A score of 4-7 denotes a usable document; 8.5 denotes professional presentation.

## Round 1 — 7.0/10

The judge inspected all six rendered pages. It found no clipping, corrupted
glyphs, missing pages, or obvious table collisions. It considered the hierarchy,
restrained colors, margins, and explicit limitations strong. The main weakness
was that comparisons concealed distinctions between methods.

The complete issue inventory and requested repairs were:

1. **Pages 2-5: curves obscured each other.** Nearly identical solid styles hid
   orange OFDG under teal gated OFDG and purple OEDG. DG disappeared in some
   smooth plots despite being listed in the legend. Shock zooms distinguished
   DG oscillations but not the filtered methods. Use consistent styles and
   sparse distinct markers; quantify coincidence with differences or ratios.
   Another physical zoom alone does not reveal intrinsically tiny differences.
2. **Page 5: cuts used different physical heights on different grids.** The
   footnote disclosed this, but the combined six-line plot invited a confounded
   refinement comparison. Evaluate at the same physical y, or separate grids
   and identify each actual y without implying a matched-cut comparison.
3. **Page 4: insufficient Shu-Osher wave detail.** The shock zoom did not examine
   the downstream oscillations motivating this experiment. Add a wave-train
   zoom and a local amplitude/detail observation, retain the shock zoom, and
   explicitly label the overview as cropped from the full domain.
4. **Pages 2, 4, 5: legends too small.** Five-entry legends in narrow columns
   were worst; repeated legends also occupied useful data regions. Use shared
   external legends and approximately 8.5-9 pt minimum final-size text.
5. **Page 2: geometry encoding was incomplete.** The six-curve straight/curved
   plot's legend named methods but not geometry. Split geometry panels or
   supply both encodings directly in the plot.
6. **Page 2: convergence interpretation needed prominence.** High filtered
   slopes could dominate a quick reading. State that they do not establish
   asymptotic order, identify nominal order, and include finest-grid errors
   or ratios alongside slopes.
7. **Page 3: false visual precision in excursions.** Tiny values such as
   4.18e-08 competed with material DG excursions. Use a stated display
   threshold and consistent significant figures; retain exact data. Identify
   the statistics immediately as sampled-polynomial excursions beyond the
   initial-data range, rather than relying on a dense later caveat.
8. **Page 4: reference sensitivity lacked an immediate consequence.** Large
   maximum differences, especially 0.463 for Shu-Osher, needed explanation
   alongside the much smaller means. Explain localized steep-front mismatch
   and limits on pointwise reference judgments. Use readable problem names.
9. **Page 5: possible color clipping was unresolved.** Determine whether the
   common map scale clips any means. If it does, state actual extrema and
   display limits with colorbar extensions; otherwise remove the vague caveat.
10. **Page 1: acronyms and notation were undefined.** Explain DG, OFDG, OEDG,
    the role of KXRCF gating, polynomial degree Pk, and cell count N. Include
    the missing Burgers grid sizes in the configuration table.
11. **Pages 1 and 6: conclusions were weaker than the available evidence.**
    Add measured takeaways on smooth error, scalar shock coincidence versus
    DG excursions, and Shu-Osher wave preservation. Tie next experiments to
    those findings. Clarify that future Riemann boundary conditions should be
    matched across comparisons; the current case already uses outflow.
12. **Page 6: dense references and vague navigation.** Separate references,
    format them consistently, add verified DOI links, and identify exact
    dataset/documentation locations. A truncated fingerprint is less useful
    to readers than an identifiable run and its manifest location.

## Revisions submitted for round 2

- Consistent method colors, solid/dashed/dotted/dash-dot lines and distinct
  sparse markers; shared external legends and larger plot text.
- Separate straight and curved smooth panels, nested markers for coincident
  errors, finest-grid error ratios, nominal orders, and prominent slope caveat.
- Filtered-method difference panels beneath scalar shock zooms; negligible
  excursions grouped below 1e-6 for display, with full values retained.
- Wider Euler panels, explicit cropped Shu-Osher overview, wave-region zoom,
  measured local density range, and retained shock/equal-DOF detail.
- Both 2D grids evaluated at y=0.740 from their saved tensor-node polynomial
  data, independently within each cell. An analytic tensor polynomial and
  interpolation at the saved nodes validate the extraction. No smoothing
  across DG interfaces or new solver runs were used.
- 2D grid resolutions shown separately; maximum gated/ungated cut differences
  stated numerically. One shared map scale includes all displayed cell means.
- Embedded text fonts, clearer definitions, completed configuration table,
  measured takeaways, evidence-linked next experiments, precise relative
  paths, and DOI links verified against the local paper title pages.

Review PDF snapshots are retained under
`measurements/study/preview/review/round-1.pdf` and `round-2.pdf`.
The underlying 103 solver runs and four references are unchanged.

## Round 2 — 8.6/10: stop condition met

The independent judge inspected all six fresh pages and verified render
freshness. It concluded that professional standard had been reached and
explicitly instructed that no third cycle was needed. It found no clipping
or layout collisions and confirmed resolution of the major first-round
issues: method visibility, scalar difference diagnostics, separated geometry
plots, wave detail, matched physical cuts, definite color coverage, definitions,
takeaways, and navigable references.

The complete remaining minor-issue inventory is retained below. These were
not blockers to 8.5, and no further PDF revisions were made after this score.

1. **Page 3, Burgers difference:** the OEDG excursion of roughly 0.24 sets the
   vertical scale, hiding the roughly 0.00049 gated-OFDG difference near zero.
   The adjacent numerical statement makes this acceptable. In an expanded
   report, give the gated difference its own inset or panel.
2. **Pages 3 and 5:** denser than pages 1 and 6, with less breathing room in
   the lower captions/tables. Still readable within six pages. A longer
   report could move detailed diagnostics to an appendix.
3. **Page 4, equal-DOF panel:** the local degree legend reuses colors that
   identify methods elsewhere. The explicit P1/P2/P3 legend and caption
   prevent material confusion. A future version could use a separate degree
   palette or direct curve labels.
4. **Throughout:** figure/panel numbers would make verbal discussion and
   later written feedback easier, although their absence is workable in a
   short brief.
5. **Page 6:** “No result here requires a giant test to interpret” is slightly
   conversational. A future neutral alternative is “The present findings can
   be investigated further with targeted refinement.”

Final outcome: **7.0 → 8.6 in two review rounds**. The approved final PDF is
`output/pdf/supervisor-preview.pdf`, identical to the round-2 snapshot.
No additional solver runs were used for this review and revision process.
