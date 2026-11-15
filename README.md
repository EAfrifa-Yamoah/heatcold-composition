# From cold to heat: the changing composition of temperature-attributable mortality in 204 countries and territories, 1990 to 2023

Secondary analysis of Global Burden of Disease Study 2023 (GBD 2023) estimates of deaths attributable
to high and low temperature. Author: Eben Afrifa-Yamoah (Edith Cowan University).

Status: first complete draft (15 November 2026). The numbers in the manuscript have been produced by
the scripts in this repository and checked by a second code path (`R/99_independent_check.R`), but
they have **not yet been reproduced by the author on his own machine**, and the references have not
yet been retrieved and read by the author. Both steps are required before submission (see
`docs/before_submission.md`).

## Layout

| Path | Content |
|---|---|
| `data/raw/` | 222 extracts of GBD 2023 estimates exactly as returned by the IHME GBD 2023 data connector on 15 November 2026 (`_harvest_log.csv` lists them). Files ending `jsonrounded` are at GBD publication precision (three or more significant figures); all others are unrounded. |
| `data/panel/` | Analysis panels built by `R/01_build.R`; `_build_report.txt` holds the build checks. Panels ending `_adj` exclude the outcome exposure to forces of nature (GBD cause 729). |
| `R/00_harvest.R` | Names the saved connector responses from their own content. Needs the session's tool-results folder and is included for transparency; it does not need to be rerun. |
| `R/01_build.R` | Builds the panels, the adjusted series and the list of shock years. |
| `R/02_analysis.R` | All analyses. Writes `output/stats.txt` (every statistic quoted in the manuscript) and `output/*.csv`. |
| `R/03_figures.R`, `R/04_tables.R` | Figures (`figures/`) and LaTeX tables (`manuscript/tables/`). |
| `R/99_independent_check.R` | Recomputes headline results from `data/raw/` with base R only and stops if any differs from the pipeline. |
| `manuscript/` | Not included in this repository. `R/04_tables.R` recreates `manuscript/tables/` locally when run; the folder is ignored by git. |
| `docs/` | Reference verification log and the list of steps the author must complete before submission. |

## Reproducing the results

R 4.3.3 with dplyr, tidyr, readr, purrr, stringr, ggplot2, patchwork, scales, sf, rnaturalearthdata,
ggrepel and jsonlite. From the repository root:

```
Rscript R/01_build.R
Rscript R/02_analysis.R          # about 20 seconds; seeds are fixed
Rscript R/99_independent_check.R
Rscript R/03_figures.R
Rscript R/04_tables.R
```

`output/stats.txt` should be identical to the committed copy apart from its time stamp.

## Things to know

* The all cause heat and cold estimates published by GBD include deaths from natural disasters
  (cause 729, exposure to forces of nature). The primary analysis removes that outcome and sets
  aside 27 country years in which it exceeds 5% of all deaths (`data/panel/shock_years.csv`).
  Section 12 of `output/stats.txt` repeats the headline results on the series as published.
* The explanation offered in the manuscript for the behaviour of the estimates in disaster years
  (a year-specific minimum risk temperature that follows the cause mix) is an inference from the
  published outputs and the methods appendix. It has not been confirmed with IHME.
* Intervals for derived indicators are approximations from published point estimates and intervals,
  not from GBD draws.
* The GBD temperature covariate is carried forward unchanged after 2020; analyses that use it stop
  in 2020.

## Use of generative AI

Claude (Anthropic) was used to retrieve the estimates through the GBD 2023 connector, to write and
run the R code, to search for and check references, and to draft the manuscript. See the declaration
in the manuscript, which must be completed by the author after verification.

## Licence

Code: MIT (see `LICENSE`). GBD estimates are the property of the Institute for Health Metrics and
Evaluation and are redistributed here for reproducibility under the IHME free-of-charge
non-commercial user agreement; check those terms before making the repository public.
