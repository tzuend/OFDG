# Report revision completed, 7 September 2026

The report was revised directly in `report/` in the codex_dev repository.
The root agent completed the editorial revision. No subagent was restarted and
no independent rating is assigned to this version. A later reviewer should
assess this compiled version rather than the earlier, partially edited draft.

## Changes

- Introduced DG through the FEM connection, independent traces, numerical flux,
  weak form and cell balance before explaining stabilization.
- Defined mapped reference spaces and successively projected physical derivatives
  before using them in curved OFDG sensor formulas.
- Qualified invariance, conservation, geometric support, cadence, and positivity
  claims. Kept the native NURBS limitation separate from polynomial curved support.
- Moved coefficient/operator details to Appendix A.
- Replaced the obsolete short development-results template with the measured
  supervisor preview: completed configurations, smooth error/rate table, scalar
  shocks, Euler shock and wave zooms, and matched physical 2D cuts.
- Checked the actual driver configuration: 2D transport velocity is (0.7,0.3);
  Lax/Shu–Osher have fixed exterior boundary data; 2D Riemann uses outflow.
- Rewrote the discussion and conclusion around supported preliminary findings.
- Added four visible, identified missing-work boxes and an appendix register.
  REF-01: reference accuracy; KX-01: spatial detection/thresholds; TIME-01:
  cadence attribution; GEO-01: broader geometric robustness.
- Added an isolated report export with input hashes. `make report-draft` reuses
  existing measurements and does not launch simulations or replace the supervisor
  brief. It also refreshes `output/pdf/thesis-revised.pdf`.

## Validation

The draft LaTeX/BibTeX build passes. The final log has no unresolved references,
overflowing boxes, unstable labels, or duplicate PDF destinations. The report
analysis agrees exactly with the retained preview analysis. All 59 PDF pages
were rendered and visually inspected using page contact sheets, with detailed
inspection of plot pages and changed layout pages. The unsigned declaration
form and the template's intentional recto/verso blank pages are retained.
`git diff --check` passes. No new numerical simulation was run for this revision.

This is a complete editorial revision of a preliminary research report, not a
claim that the explicitly listed research gaps have been resolved. The final
publication guard remains in place.

Delivered PDF SHA-256: `96c041e56432dc085fc928a6b4f1dd93e7e68e5031bd78f1132ab96856ae22d8`

## Curvilinear mathematics addition

Section 4.2 (printed pages 20–23) now derives physical weighted mass and
nested projection matrices, the weak projected derivative solve, recursive
higher derivatives and their projection defect, and curved-face metrics and
normalization. A one-dimensional mapped example explains why the derivative
leaves the represented space. Mean preservation and recovery of the affine
construction are stated with the quadrature qualifications. No simulations
or subagent review were run. The expanded PDF was rendered and inspected;
the contents spacing was adjusted to avoid an isolated bibliography entry.
