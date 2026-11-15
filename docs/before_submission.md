# Steps the author must complete before submission

1. Rerun `R/01_build.R` to `R/04_tables.R` and `R/99_independent_check.R` on your own machine and
   confirm that `output/stats.txt` is unchanged. Check every number in `manuscript/main.tex` against
   `output/stats.txt` (sections are numbered to match the order of the Results).
2. Spot check the extracts in `data/raw/` against the GBD Results Tool
   (https://vizhub.healthdata.org/gbd-results), in particular: global heat and cold deaths in 1990,
   2010 and 2023; Haiti 2010 (high temperature, all causes; exposure to forces of nature;
   cardiovascular diseases); Libya 2023. The extracts were retrieved through the IHME GBD 2023
   connector and have not been compared with the Results Tool.
3. Retrieve and read every reference (`docs/reference_verification.md` records what was and was
   not confirmed during drafting).
4. Run and document the literature search for the Research in context panel. The panel currently
   says that this has not been done.
5. Consider writing to IHME about the disaster-year behaviour before submission. The mechanism
   proposed in the manuscript is an inference and IHME can confirm or correct it.
6. Decide the threshold for shock years (5% of all deaths is a judgement) and, if changed, set
   `SHOCK_SHARE` in `R/01_build.R` and rerun.
7. Complete the AI-use declaration, the contributors statement and the repository URL, all marked
   in red in the PDF, and revise the prose towards your own voice.
8. Check the IHME data use terms before making `data/raw/` public.
