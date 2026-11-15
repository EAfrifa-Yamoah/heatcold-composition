# 04_tables.R
# LaTeX tables for the manuscript and supplement, written to manuscript/tables/ from output/*.csv
# and data/panel/*.csv. Run after R/02_analysis.R.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)
ro <- function(f) read_csv(file.path("output", f), show_col_types = FALSE)
n0 <- function(x) formatC(round(x), format = "d", big.mark = ","); f1 <- function(x, d = 1) formatC(x, format = "f", digits = d)
tex <- function(x) gsub("&", "\\\\&", gsub("ü", "\\\\\"u", gsub("ô", "\\\\^o", x)))
wr <- function(lines, f) writeLines(lines, file.path("manuscript/tables", f))
reg_order <- c("Global", "High income", "Central and eastern Europe, central Asia", "Latin America and Caribbean", "North Africa and Middle East", "South Asia", "Southeast and east Asia, Oceania", "Sub-Saharan Africa")

# Table 1: levels by super region
t1 <- ro("table1_regions.csv") %>% mutate(region = factor(region, reg_order)) %>% arrange(region, year)
rows <- t1 %>% transmute(a = ifelse(year == 1990, as.character(region), ""), year, heat = sprintf("%s (%s to %s)", n0(h_n), n0(h_n_lo), n0(h_n_hi)), cold = sprintf("%s (%s to %s)", n0(c_n), n0(c_n_lo), n0(c_n_hi)),
                         share = f1(share_n), hr = f1(h_asr, 2), cr = f1(c_asr, 2), hf = f1(100 * h_asr / m_asr, 2), cf = f1(100 * c_asr / m_asr, 2))
wr(c("\\begin{tabular}{llrrrrrrr}", "\\toprule", " & & \\multicolumn{2}{c}{Attributable deaths (95\\% interval)} & Heat share & \\multicolumn{2}{c}{Rate per 100,000} & \\multicolumn{2}{c}{Fraction (\\%)} \\\\",
     "\\cmidrule(lr){3-4}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}", " & Year & Heat & Cold & (\\%) & Heat & Cold & Heat & Cold \\\\", "\\midrule",
     apply(rows, 1, function(r) paste(paste(r, collapse = " & "), "\\\\", ifelse(r["year"] == "2023" , "\\addlinespace[2pt]", ""))), "\\bottomrule", "\\end{tabular}"), "table1_regions.tex")

# Table 2: decomposition of the change in counts
d <- ro("decomposition_counts.csv") %>% mutate(region = factor(region, reg_order), risk = ifelse(risk == "h", "Heat", "Cold")) %>% arrange(desc(risk), region)
sg <- function(x) sprintf("%+.1f", x)
rows <- d %>% transmute(a = as.character(region), n90 = n0(n90), n23 = n0(n23), pop = sg(pop_pct), age = sg(age_pct), mort = sg(mort_pct), paf = sg(paf_pct), net = sg(net_pct))
body <- c("\\multicolumn{8}{l}{\\textit{Heat}} \\\\", apply(rows[d$risk == "Heat", ], 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), "\\addlinespace[3pt]", "\\multicolumn{8}{l}{\\textit{Cold}} \\\\", apply(rows[d$risk == "Cold", ], 1, function(r) paste(paste(r, collapse = " & "), "\\\\")))
wr(c("\\begin{tabular}{lrrrrrrr}", "\\toprule", " & \\multicolumn{2}{c}{Attributable deaths} & \\multicolumn{5}{c}{Contribution to the change (\\% of the 1990 count)} \\\\", "\\cmidrule(lr){2-3}\\cmidrule(lr){4-8}",
     " & 1990 & 2023 & Population & Age & All cause & Attributable & Net \\\\", " & & & growth & structure & mortality & fraction & change \\\\", "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "table2_decomposition.tex")

# Table 3: shock years
sy <- read_csv("data/panel/shock_years.csv", show_col_types = FALSE) %>% filter(d729_n >= 1000) %>% left_join(ro("shock_cases.csv") %>% select(location_id, year, cvd_heat, cvd_heat_adjacent, cvd_cold, cvd_cold_adjacent), by = c("location_id", "year")) %>% arrange(desc(d729_share))
cvf <- function(a, b) ifelse(is.na(a), "..", sprintf("%s (%s)", n0(a), n0(b)))
rows <- sy %>% transmute(a = tex(sub(" \\(.*", "", location_name)), year, d = n0(d729_n), sh = f1(d729_share), hp = n0(heat_published), h7 = n0(heat_729), pf = f1(heat_729_paf), cp = n0(cold_published), c7 = n0(cold_729), cvh = cvf(cvd_heat, cvd_heat_adjacent), cvc = cvf(cvd_cold, cvd_cold_adjacent))
wr(c("\\begin{tabular}{llrrrrrrrrr}", "\\toprule", " & & \\multicolumn{2}{c}{Forces of nature deaths} & \\multicolumn{3}{c}{Heat attributable deaths} & \\multicolumn{2}{c}{Cold attributable deaths} & \\multicolumn{2}{c}{Cardiovascular deaths attributed to} \\\\",
     "\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}", "Country & Year & Number & \\% of all & All causes, & Forces of & Fraction & All causes, & Forces of & Heat & Cold \\\\",
     " & & & deaths & published & nature & (\\%) & published & nature & (adjacent years) & (adjacent years) \\\\", "\\midrule", apply(rows, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), "\\bottomrule", "\\end{tabular}"), "table3_shock_years.tex")

# Supplementary: causes
gw <- ro("global_causes.csv"); sh <- ro("shapley_share_causes.csv") %>% select(side, cause, points) %>% pivot_wider(names_from = side, values_from = points)
rows <- gw %>% left_join(sh, by = "cause") %>% transmute(a = cause, w0 = f1(100 * w_1990), w1 = f1(100 * w_2023), h0 = n0(h_1990), h1 = n0(h_2023), c0 = n0(c_1990), c1 = n0(c_2023), ph0 = f1(100 * ph_1990, 2), ph1 = f1(100 * ph_2023, 2), pc0 = f1(100 * pc_1990, 2), pc1 = f1(100 * pc_2023, 2), sh_h = sg(heat), sh_c = sg(cold))
wr(c("\\begin{tabular}{lrrrrrrrrrrrr}", "\\toprule", " & \\multicolumn{2}{c}{Share of all deaths (\\%)} & \\multicolumn{2}{c}{Heat deaths} & \\multicolumn{2}{c}{Cold deaths} & \\multicolumn{2}{c}{Heat fraction (\\%)} & \\multicolumn{2}{c}{Cold fraction (\\%)} & \\multicolumn{2}{c}{Shapley (points)} \\\\",
     "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\\cmidrule(lr){12-13}", "Cause group & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 & Heat & Cold \\\\", "\\midrule",
     apply(rows, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), "\\bottomrule", "\\end{tabular}"), "tableS_causes.tex")

# Supplementary: all shock years
sa <- read_csv("data/panel/shock_years.csv", show_col_types = FALSE) %>% arrange(desc(d729_share))
rows <- sa %>% transmute(a = tex(sub(" \\(.*", "", location_name)), year, d = n0(d729_n), sh = f1(d729_share), hp = f1(heat_published), h7 = f1(heat_729), pf = f1(heat_729_paf), cp = f1(cold_published), c7 = f1(cold_729))
wr(c("\\begin{tabular}{llrrrrrrr}", "\\toprule", "Country & Year & Forces of nature & \\% of all & Heat, all causes & Heat, forces & Fraction & Cold, all causes & Cold, forces \\\\", " & & deaths & deaths & (published) & of nature & (\\%) & (published) & of nature \\\\", "\\midrule",
     apply(rows, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")), "\\bottomrule", "\\end{tabular}"), "tableS_shock_all.tex")

# Supplementary: country table (longtable)
ct <- ro("country_table.csv") %>% mutate(region = factor(region, reg_order[-1])) %>% arrange(region, location_name)
cl <- function(x) ifelse(is.finite(x) & x >= -50 & x <= 150, f1(x), "n/a")
rows <- ct %>% transmute(a = paste0(tex(location_name), ifelse(y0 != 1990 | y1 != 2023, "$^{\\dagger}$", "")), h0 = n0(heat_n_1990), h1 = n0(heat_n_2023), c0 = n0(cold_n_1990), c1 = n0(cold_n_2023), hr0 = f1(heat_rate_1990, 2), hr1 = f1(heat_rate_2023, 2), cr0 = f1(cold_rate_1990, 2), cr1 = f1(cold_rate_2023, 2),
                         s0 = cl(share90), s1 = ifelse(is.finite(share23) & share23 >= -50 & share23 <= 150, sprintf("%s (%s to %s)", f1(share23), f1(pmax(lo23, -99)), f1(pmin(hi23, 199))), "n/a"), pr = f1(p_rise_shared, 2))
lines <- c("\\begin{longtable}{p{4.3cm}rrrrrrrrrlr}", "\\toprule", " & \\multicolumn{2}{c}{Heat deaths} & \\multicolumn{2}{c}{Cold deaths} & \\multicolumn{2}{c}{Heat rate} & \\multicolumn{2}{c}{Cold rate} & \\multicolumn{2}{c}{Heat share (\\%)} & P(rise) \\\\",
           "Country or territory & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 & 1990 & 2023 (95\\% interval) & \\\\", "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (r in levels(ct$region)) { lines <- c(lines, sprintf("\\multicolumn{12}{l}{\\textit{%s}} \\\\", r), apply(rows[ct$region == r, ], 1, function(x) paste(paste(x, collapse = " & "), "\\\\")), "\\addlinespace[3pt]") }
wr(c(lines, "\\end{longtable}"), "tableS_countries.tex")
cat("tables written:", paste(list.files("manuscript/tables"), collapse = ", "), "\n")
