# Suggested references for supervisor review

The active thesis bibliography contains only the three research papers and the Hesthaven–Warburton textbook excerpt currently in `papers/`. All other existing bibliography entries were moved to `suggested.bib`; no additional PDFs were found to move, and no papers were downloaded. These are candidates, not an instruction to cite everything. Metadata inherited from the old bibliography should be checked against each paper before inclusion.

## Most directly relevant

| Entry | Why it may belong in the thesis |
|---|---|
| Krivodonova et al. (2004), `KrivodonovaEtAl2004` | Original KXRCF source. Highest priority: needed to establish attribution and distinguish the implemented radius, component pooling, and inflow choices from the original method. |
| Zhang & Shu (2010), `ZhangShu2010` | Mean-preserving positivity scaling for Euler DG. Needed to explain the origin and assumptions of the separate admissibility control, without transferring a theorem to an unverified configuration. |
| Anderson et al. (2021), `AndersonEtAl2021` | Scholarly attribution for MFEM and context for the finite-element implementation. Useful even if the thesis is not a software study. |
| Cockburn & Shu (1998), `CockburnShu1998` | Foundational RKDG treatment for multidimensional systems. Helpful if the Methods chapter needs more detailed support than the supplied textbook excerpt and OFDG papers. |

## Conditional on the final questions and experiments

| Entry | Reason and condition |
|---|---|
| Qiu & Shu (2005), `QiuShu2005` | Useful if comparing troubled-cell indicators or discussing detector sensitivity. |
| Fu & Shu (2017), `FuShu2017` | Alternative detector formulation; useful for a focused comparison, not required merely to implement KXRCF. |
| Vuik & Ryan (2014), `VuikRyan2014` | Multiwavelet detection, relevant only if the scope broadens to detector alternatives. |
| Vuik & Ryan (2016), `VuikRyan2016` | Automatic parameter selection; relevant if threshold sensitivity becomes a research question. |
| Jiang & Shu (1996), `JiangShu1996` | WENO reference-solver reconstruction, if independent WENO solutions support the final Results chapter. |
| Gottlieb & Shu (1998), `GottliebShu1998` | TVD/SSP time-stepping background for the reference solver or a cadence study. It is not a citation for classical RK4 being SSP. |
| Sod (1978), `Sod1978` | Reference-solver shock-tube verification, if that verification is documented in the thesis. |
| Shu & Osher (1989), `ShuOsher1989` | Original context for a shock–entropy-wave experiment if retained. |
| Lax & Liu (1998), `LaxLiu1998` | Source context for a final two-dimensional Riemann experiment. |
| Woodward & Colella (1984), `WoodwardColella1984` | Relevant if strong-shock/blast-wave problems enter the final study. Low priority for the current limited scope. |

## Bibliography corrections made now

- Removed an unrelated radiative-transfer abstract accidentally stored in the 2021 OFDG entry.
- Removed placeholder abstract metadata from the OEDG entry.
- Used the supplied OEDG PDF's journal information: volume 94, issue 353, pages 1147–1198, 2025, with a note that it appeared online in 2024. The existing citation key is retained for stable source references.
- Added the Hesthaven–Warburton textbook that is already present as a chapter excerpt, rather than introducing another external DG reference.
