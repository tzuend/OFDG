# Thesis revision — 9 September 2026

The thesis now follows Introduction, Mathematical Background, Methods, Results, Discussion, Conclusion, with a provisional abstract and a coefficient-space appendix. Mathematical Background is a separate input file so it can be relocated after supervisor feedback.

The Methods chapter contains the adapted OFDG definition, KXRCF, physical curvilinear projections and recursive projected derivatives, affine OEDG, numerical controls, and implementation assumptions. Empirical preview conclusions, timing statistics, and validation-result tables have been removed from the compiled thesis. Results, Discussion, and Conclusion each contain a visible requirements placeholder. The abstract is explicitly provisional. Source gaps for original KXRCF and positivity theory are also marked.

Only the three research papers and the Hesthaven–Warburton textbook excerpt in papers/ are cited. Other entries are retained in suggested-references/suggested.bib with selection notes in suggested-references/README.md. No external PDFs were downloaded or found to move. Bibliography metadata corrections are documented there.

Build from the repository root with `make report-draft`. This compiles the thesis without regenerating or importing preliminary experiment assets. The separate supervisor-preview PDF and raw numerical measurements are preserved. The editable source is report/thesis.tex and its included files; the output is output/pdf/thesis-revised.pdf.

The previous draft sources and PDF are preserved in tmp/thesis-before-methods-revision-20260909/. The final draft contains 42 pages including front matter, intentional recto blanks, and the declaration form. Every page was visually inspected. Citation, placeholder, and clean-build checks are recorded in tmp/thesis-methods-qa/delivery/audit.json. No numerical solvers or reviewer subagents were run for this revision.
