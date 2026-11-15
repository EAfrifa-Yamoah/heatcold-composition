# 02_analysis.R
# "From cold to heat: the changing composition of temperature-attributable mortality in 204
# countries, 1990 to 2023". Every statistic quoted in the manuscript is written to
# output/stats.txt; tables for figures go to output/*.csv. Run from the repository root after
# R/01_build.R:  Rscript R/02_analysis.R
#
# Primary series: heat and cold attributable deaths excluding the outcome "exposure to forces of
# nature" (cause 729), the "adjusted" panels of R/01_build.R. Shock years (cause 729 at least 5% of
# all deaths in a country year) are excluded from country level analyses; at the two ends of the
# period the nearest year that is not a shock year is used (columns y0 and y1 of endpoints.csv).
# Section 12 repeats the headline results on the series as published.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr)})
set.seed(20261115); dir.create("output", showWarnings = FALSE)
rp <- function(f) read_csv(file.path("data/panel", f), show_col_types = FALSE)
series_pub <- rp("series.csv"); series <- rp("series_adj.csv"); agg_pub <- rp("agg_series.csv"); agg <- rp("agg_series_adj.csv")
byage <- rp("byage_adj.csv"); causes <- rp("causes_adj.csv"); parent <- rp("parent_adj.csv"); yll <- rp("yll_adj.csv")
temp <- rp("temperature.csv"); ctry <- rp("countries.csv"); ep <- rp("endpoints.csv"); shock <- rp("shock_years.csv"); dis <- rp("disaster.csv")
S <- series %>% filter(sex == "combined"); Sp <- series_pub %>% filter(sex == "combined")
NDRAW <- 4000
sr_short <- c("Central Europe, Eastern Europe, and Central Asia" = "Central and eastern Europe, central Asia", "High-income" = "High income",
              "Latin America and Caribbean" = "Latin America and Caribbean", "North Africa and Middle East" = "North Africa and Middle East",
              "South Asia" = "South Asia", "Southeast Asia, East Asia, and Oceania" = "Southeast and east Asia, Oceania", "Sub-Saharan Africa" = "Sub-Saharan Africa", "Global" = "Global")

# ---- helpers ---------------------------------------------------------------------------------
# exact Shapley values by subset enumeration; x0, x1 are lists of factors (scalars or vectors)
shapley <- function(f, x0, x1) {
  k <- length(x0); v <- numeric(2^k)
  for (mask in 0:(2^k - 1)) { on <- bitwAnd(mask, 2^(0:(k - 1))) > 0; x <- x0; x[on] <- x1[on]; v[mask + 1] <- f(x) }
  phi <- numeric(k)
  for (i in 1:k) for (mask in 0:(2^k - 1)) if (!bitwAnd(mask, 2^(i - 1))) {
    s <- sum(bitwAnd(mask, 2^(0:(k - 1))) > 0)
    phi[i] <- phi[i] + factorial(s) * factorial(k - s - 1) / factorial(k) * (v[mask + 2^(i - 1) + 1] - v[mask + 1])
  }
  stopifnot(abs(sum(phi) - (v[2^k] - v[1])) < 1e-8 * max(1, abs(v[2^k]))); setNames(phi, names(x0))
}
# two piece normal and split log normal draws reproducing a point estimate and 95% interval
tpn  <- function(z, m, lo, hi) m + z * ifelse(z < 0, (m - lo) / 1.96, (hi - m) / 1.96)
sln  <- function(z, m, lo, hi) m * exp(z * ifelse(z < 0, (log(m) - log(lo)) / 1.96, (log(hi) - log(m)) / 1.96))
rcop <- function(n, rho) { z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n); cbind(z1, z2) }
qi   <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
fmt  <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
pct  <- function(a, b) 100 * (b / a - 1)
gini <- function(x) { x <- sort(x); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }
# country values at the two ends of the period (y0, y1), with suffixes _1990 and _2023
ends <- function(d, vars) {
  a <- d %>% inner_join(ep %>% select(location_id, y0), by = "location_id") %>% filter(year == y0) %>% select(location_id, all_of(vars)) %>% rename_with(~ paste0(.x, "_1990"), all_of(vars))
  b <- d %>% inner_join(ep %>% select(location_id, y1), by = "location_id") %>% filter(year == y1) %>% select(location_id, all_of(vars)) %>% rename_with(~ paste0(.x, "_2023"), all_of(vars))
  ctry %>% inner_join(a, by = "location_id") %>% inner_join(b, by = "location_id")
}

sink("output/stats.txt")
cat("From cold to heat: statistics\nGenerated", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Countries:", n_distinct(S$location_id), " years:", paste(range(S$year), collapse = " to "), "\n")

# ==== 0. The disaster outcome and shock years ======================================================
cat("\n==== 0. Exposure to forces of nature (cause 729) inside the published heat and cold estimates ====\n")
gp <- agg_pub %>% filter(location_id == 1) %>% select(year, hp_n = h_n, cp_n = c_n) %>% inner_join(agg %>% filter(location_id == 1) %>% select(year, h729 = h729_n_removed, c729 = c729_n_removed, d729 = d729_n_removed), by = "year") %>%
  mutate(h_share = 100 * h729 / hp_n, c_share = 100 * c729 / cp_n)
cat(sprintf("Global heat attributable deaths from cause 729, 1990 to 2023: total %s; median year %s; five largest years: %s\n", fmt(sum(gp$h729), 0), fmt(median(gp$h729), 0),
            paste(sprintf("%d %s (%.1f%% of published heat deaths)", gp$year[order(-gp$h729)][1:5], fmt(sort(gp$h729, TRUE)[1:5], 0), gp$h_share[order(-gp$h729)][1:5]), collapse = "; ")))
cat(sprintf("Global cold attributable deaths from cause 729: negative in %d of 34 years; total %s; largest in absolute value: %s\n", sum(gp$c729 < 0), fmt(sum(gp$c729), 0),
            paste(sprintf("%d %s", gp$year[order(gp$c729)][1:3], fmt(sort(gp$c729)[1:3], 0)), collapse = "; ")))
cat(sprintf("1990: heat component %s (%.1f%% of published), cold component %s; 2023: heat %s (%.1f%%), cold %s\n", fmt(gp$h729[1], 0), gp$h_share[1], fmt(gp$c729[1], 0), fmt(gp$h729[34], 0), gp$h_share[34], fmt(gp$c729[34], 0)))
D <- Sp %>% select(location_id, location_name, year, h_n, c_n, d_n) %>% inner_join(dis %>% filter(sex == "combined"), by = c("location_id", "year")) %>% mutate(share = 100 * d729_n / d_n)
cat(sprintf("Shock years (cause 729 at least 5%% of all deaths): %d country years in %d countries; of the 6,936 country years, %d have a share of at least 1%%\n", nrow(shock), n_distinct(shock$location_id), sum(D$share >= 1)))
cat(sprintf("Country years in which the heat attributable deaths from cause 729 exceed the published heat total for all causes (deaths above 1): %d (%s)\n", sum(D$h729_n > D$h_n & D$h729_n > 1),
            paste(sprintf("%s %d", D$location_name[D$h729_n > D$h_n & D$h729_n > 1], D$year[D$h729_n > D$h_n & D$h729_n > 1]), collapse = "; ")))
pafq <- D %>% filter(d729_n >= 10) %>% mutate(grp = ifelse(share >= 5, "shock", "other"))
cat(sprintf("Heat attributable fraction of cause 729 deaths (country years with at least 10 such deaths): median %.1f%% (IQR %.1f to %.1f) in %d non-shock country years; %.1f%% (IQR %.1f to %.1f) in %d shock country years\n",
            median(pafq$h729_paf[pafq$grp == "other"]), quantile(pafq$h729_paf[pafq$grp == "other"], .25), quantile(pafq$h729_paf[pafq$grp == "other"], .75), sum(pafq$grp == "other"),
            median(pafq$h729_paf[pafq$grp == "shock"]), quantile(pafq$h729_paf[pafq$grp == "shock"], .25), quantile(pafq$h729_paf[pafq$grp == "shock"], .75), sum(pafq$grp == "shock")))
# cardiovascular attribution around selected events (12 countries extracted year by year)
cv <- read_csv(list.files("data/raw", "multi_location_trend__r337__c491__deaths__all__1990to2023", full.names = TRUE), show_col_types = FALSE) %>% transmute(location_id, location_name, year, cvd_heat = count_combined, cvd_heat_paf = percent_combined) %>%
  inner_join(read_csv(list.files("data/raw", "multi_location_trend__r338__c491__deaths__all__1990to2023", full.names = TRUE), show_col_types = FALSE) %>% transmute(location_id, year, cvd_cold = count_combined, cvd_cold_paf = percent_combined), by = c("location_id", "year"))
ev <- shock %>% filter(location_id %in% cv$location_id) %>% select(location_id, year, d729_share, heat_729_paf)
cases <- map_dfr(seq_len(nrow(ev)), function(i) { x <- cv %>% filter(location_id == ev$location_id[i]); e <- x %>% filter(year == ev$year[i]); nb <- x %>% filter(year %in% (ev$year[i] + c(-1, 1)))
  tibble(location_id = e$location_id, location_name = e$location_name, year = e$year, share_729 = ev$d729_share[i], heat_729_paf = ev$heat_729_paf[i], cvd_heat = e$cvd_heat, cvd_heat_adjacent = mean(nb$cvd_heat), cvd_heat_paf = e$cvd_heat_paf, cvd_heat_paf_adjacent = mean(nb$cvd_heat_paf),
         cvd_cold = e$cvd_cold, cvd_cold_adjacent = mean(nb$cvd_cold), cvd_cold_paf = e$cvd_cold_paf, cvd_cold_paf_adjacent = mean(nb$cvd_cold_paf)) }) %>% arrange(desc(share_729))
write_csv(cases, "output/shock_cases.csv")
cat("Cardiovascular deaths attributed to heat and to cold in shock years against the mean of the adjacent years (rounded source values):\n")
print(as.data.frame(cases %>% transmute(location = substr(location_name, 1, 26), year, share_729 = round(share_729, 1), fon_heat_paf = round(heat_729_paf, 1), cvd_heat = round(cvd_heat), adj = round(cvd_heat_adjacent), cvd_heat_paf = round(cvd_heat_paf, 2), adj_paf = round(cvd_heat_paf_adjacent, 2),
                                      cvd_cold = round(cvd_cold), adj_cold = round(cvd_cold_adjacent))), row.names = FALSE)
cat(sprintf("Of %d shock years with cardiovascular series: heat attributable cardiovascular deaths negative in %d; cold attributable cardiovascular deaths below half of the adjacent mean in %d\n", nrow(cases), sum(cases$cvd_heat < 0), sum(cases$cvd_cold < 0.5 * cases$cvd_cold_adjacent)))
write_csv(gp, "output/global_disaster_component.csv")

# ==== 1. Global and super region headline ========================================================
cat("\n==== 1. Global and super region levels, 1990 and 2023 (adjusted series; published in brackets where shown) ====\n")
A <- agg %>% mutate(share_n = 100 * h_n / (h_n + c_n), share_asr = 100 * h_asr / (h_asr + c_asr), region = sr_short[location_name])
Ap <- agg_pub %>% mutate(share_n = 100 * h_n / (h_n + c_n))
g <- A %>% filter(location_id == 1); g90 <- g %>% filter(year == 1990); g23 <- g %>% filter(year == 2023)
p90 <- Ap %>% filter(location_id == 1, year == 1990); p23 <- Ap %>% filter(location_id == 1, year == 2023)
cat(sprintf("Global heat deaths: %s (%s to %s) in 1990, %s (%s to %s) in 2023, change %+.1f%% [published %s and %s, %+.1f%%]\n", fmt(g90$h_n, 0), fmt(g90$h_n_lo, 0), fmt(g90$h_n_hi, 0), fmt(g23$h_n, 0), fmt(g23$h_n_lo, 0), fmt(g23$h_n_hi, 0), pct(g90$h_n, g23$h_n), fmt(p90$h_n, 0), fmt(p23$h_n, 0), pct(p90$h_n, p23$h_n)))
cat(sprintf("Global cold deaths: %s (%s to %s) in 1990, %s (%s to %s) in 2023, change %+.1f%% [published %s and %s, %+.1f%%]\n", fmt(g90$c_n, 0), fmt(g90$c_n_lo, 0), fmt(g90$c_n_hi, 0), fmt(g23$c_n, 0), fmt(g23$c_n_lo, 0), fmt(g23$c_n_hi, 0), pct(g90$c_n, g23$c_n), fmt(p90$c_n, 0), fmt(p23$c_n, 0), pct(p90$c_n, p23$c_n)))
cat(sprintf("Global non-optimal temperature deaths: %s in 1990, %s in 2023 (%+.1f%%); heat + cold over parent: %.4f (1990), %.4f (2023)\n", fmt(g90$p_n, 0), fmt(g23$p_n, 0), pct(g90$p_n, g23$p_n), (g90$h_n + g90$c_n) / g90$p_n, (g23$h_n + g23$c_n) / g23$p_n))
cat(sprintf("Global deaths from all causes other than cause 729: %s to %s (%+.1f%%)\n", fmt(g90$d_n, 0), fmt(g23$d_n, 0), pct(g90$d_n, g23$d_n)))
cat(sprintf("Age standardised rate per 100,000: heat %.2f to %.2f (%+.1f%%); cold %.2f to %.2f (%+.1f%%); all cause %.1f to %.1f (%+.1f%%)\n",
            g90$h_asr, g23$h_asr, pct(g90$h_asr, g23$h_asr), g90$c_asr, g23$c_asr, pct(g90$c_asr, g23$c_asr), g90$m_asr, g23$m_asr, pct(g90$m_asr, g23$m_asr)))
cat(sprintf("Attributable fraction of all deaths, crude (%%): heat %.3f (%.3f to %.3f) to %.3f (%.3f to %.3f), relative %+.1f%%; cold %.3f (%.3f to %.3f) to %.3f (%.3f to %.3f), relative %+.1f%%\n",
            g90$h_paf, g90$h_paf_lo, g90$h_paf_hi, g23$h_paf, g23$h_paf_lo, g23$h_paf_hi, pct(g90$h_paf, g23$h_paf), g90$c_paf, g90$c_paf_lo, g90$c_paf_hi, g23$c_paf, g23$c_paf_lo, g23$c_paf_hi, pct(g90$c_paf, g23$c_paf)))
cat(sprintf("Attributable fraction, age standardised (%%): heat %.3f to %.3f (relative %+.1f%%); cold %.3f to %.3f (relative %+.1f%%)\n", 100 * g90$h_asr / g90$m_asr, 100 * g23$h_asr / g23$m_asr, pct(g90$h_asr / g90$m_asr, g23$h_asr / g23$m_asr),
            100 * g90$c_asr / g90$m_asr, 100 * g23$c_asr / g23$m_asr, pct(g90$c_asr / g90$m_asr, g23$c_asr / g23$m_asr)))
cat(sprintf("Heat share of temperature attributable deaths: %.1f%% to %.1f%% (counts) [published %.1f%% to %.1f%%]; %.1f%% to %.1f%% (age standardised rates); cold to heat ratio %.2f to %.2f\n",
            g90$share_n, g23$share_n, p90$share_n, p23$share_n, g90$share_asr, g23$share_asr, g90$c_n / g90$h_n, g23$c_n / g23$h_n))
gy <- g %>% select(year, h_n, c_n, share_n, h_asr, c_asr, h_paf, c_paf, m_asr, d_n) %>% inner_join(Ap %>% filter(location_id == 1) %>% select(year, share_pub = share_n, h_pub = h_n, c_pub = c_n), by = "year")
cat(sprintf("Global heat share by year (adjusted): highest before 2023 in %d (%.1f%%); 1998 %.1f%%; 2010 %.1f%% [published %.1f%%]; 2016 %.1f%%; 2019 %.1f%%; 2020 %.1f%%; 2021 %.1f%%; 2022 %.1f%%\n", gy$year[which.max(gy$share_n[gy$year < 2023])], max(gy$share_n[gy$year < 2023]),
            gy$share_n[gy$year == 1998], gy$share_n[gy$year == 2010], gy$share_pub[gy$year == 2010], gy$share_n[gy$year == 2016], gy$share_n[gy$year == 2019], gy$share_n[gy$year == 2020], gy$share_n[gy$year == 2021], gy$share_n[gy$year == 2022]))
lt <- lm(share_n ~ year, gy %>% filter(year <= 2019)); cat(sprintf("Linear trend of the global heat share, 1990 to 2019: %.2f percentage points per decade (SE %.2f)\n", 10 * coef(lt)[2], 10 * summary(lt)$coefficients[2, 2]))
cat(sprintf("Global heat deaths, peak years (adjusted): %s\n", paste(sprintf("%d %s", gy$year[order(-gy$h_n)][1:6], fmt(sort(gy$h_n, TRUE)[1:6], 0)), collapse = "; ")))
tab1 <- A %>% filter(year %in% c(1990, 2023)) %>% select(region, location_id, year, h_n, h_n_lo, h_n_hi, c_n, c_n_lo, c_n_hi, h_asr, c_asr, h_paf, c_paf, share_n, share_asr, d_n, m_asr) %>%
  left_join(Ap %>% select(location_id, year, h_pub = h_n, c_pub = c_n, share_pub = share_n), by = c("location_id", "year"))
write_csv(tab1, "output/table1_regions.csv")
cat("\nSuper regions (adjusted): heat deaths, cold deaths, heat share of counts %, crude heat and cold fractions %, age standardised rates:\n")
print(as.data.frame(tab1 %>% transmute(region = substr(region, 1, 30), year, heat = round(h_n), cold = round(c_n), share = round(share_n, 1), share_pub = round(share_pub, 1), h_paf = round(h_paf, 3), c_paf = round(c_paf, 3), h_asr = round(h_asr, 2), c_asr = round(c_asr, 2)) %>% arrange(region, year)), row.names = FALSE)
write_csv(A %>% left_join(Ap %>% select(location_id, year, share_pub = share_n, h_pub = h_n, c_pub = c_n), by = c("location_id", "year")), "output/agg_shares.csv")

# ==== 2. Demographic decomposition of counts (Global and super regions) =========================
cat("\n==== 2. Decomposition of the change in attributable deaths, 1990 to 2023 ====\n")
cat("C = N * sum_a s_a m_a p_a: population size N, age structure s, age specific all cause death rate m, age specific attributable fraction p; Shapley values over the four factors\n")
dec <- map_dfr(unique(byage$location_id), function(L) map_dfr(c("h", "c"), function(r) {
  b0 <- byage %>% filter(location_id == L, year == 1990) %>% arrange(age_group_id); b1 <- byage %>% filter(location_id == L, year == 2023) %>% arrange(age_group_id)
  stopifnot(identical(b0$age_group_id, b1$age_group_id))
  mk <- function(b) list(N = sum(b$pop), s = b$pop / sum(b$pop), m = b$d_n / b$pop, p = b[[paste0(r, "_n")]] / b$d_n)
  f <- function(x) x$N * sum(x$s * x$m * x$p)
  x0 <- mk(b0); x1 <- mk(b1); phi <- shapley(f, x0, x1)
  tibble(location_id = L, location_name = b0$location_name[1], risk = r, n90 = f(x0), n23 = f(x1), pop = phi["N"], age = phi["s"], mort = phi["m"], paf = phi["p"])
})) %>% mutate(region = sr_short[location_name], across(c(pop, age, mort, paf), ~ 100 * .x / n90, .names = "{.col}_pct"), net_pct = 100 * (n23 / n90 - 1))
write_csv(dec, "output/decomposition_counts.csv")
chk <- dec %>% left_join(A %>% filter(year == 2023) %>% select(location_id, h_n, c_n), by = "location_id") %>% mutate(rel = n23 / ifelse(risk == "h", h_n, c_n) - 1)
cat(sprintf("Check: sum of age specific attributable deaths against the all age count, max absolute relative difference %.2e\n", max(abs(chk$rel))))
cat("Contributions as % of the 1990 count (population growth, age structure, all cause mortality, attributable fraction, net):\n")
print(as.data.frame(dec %>% transmute(region = substr(region, 1, 30), risk = ifelse(risk == "h", "heat", "cold"), n90 = round(n90), n23 = round(n23), pop = round(pop_pct, 1), age = round(age_pct, 1), mort = round(mort_pct, 1), paf = round(paf_pct, 1), net = round(net_pct, 1)) %>% arrange(risk, region)), row.names = FALSE)
ga <- byage %>% filter(location_id == 1); old <- c(19, 20, 30, 31, 32, 235); u5 <- c(2, 3, 388, 389, 238, 34)
agesh <- ga %>% group_by(year) %>% summarise(pop70 = 100 * sum(pop[age_group_id %in% old]) / sum(pop), popn = sum(pop),
           h_u5 = 100 * sum(h_n[age_group_id %in% u5]) / sum(h_n), h_1549 = 100 * sum(h_n[age_group_id %in% 8:14]) / sum(h_n), h_70 = 100 * sum(h_n[age_group_id %in% old]) / sum(h_n),
           c_u5 = 100 * sum(c_n[age_group_id %in% u5]) / sum(c_n), c_1549 = 100 * sum(c_n[age_group_id %in% 8:14]) / sum(c_n), c_70 = 100 * sum(c_n[age_group_id %in% old]) / sum(c_n),
           c_neg = sum(c_n[c_n < 0]), c_pos = sum(c_n[c_n > 0]))
cat("\nGlobal age pattern: % of population aged 70+, population, % of heat and of cold deaths under 5, 15 to 49, 70+; sum of negative and of positive age specific cold deaths:\n"); print(as.data.frame(agesh %>% mutate(across(-c(year, popn, c_neg, c_pos), ~ round(.x, 1)), across(c(popn, c_neg, c_pos), round))), row.names = FALSE)
write_csv(ga %>% select(year, age_group_id, age_group, h_n, c_n, h_paf, c_paf, d_n, pop), "output/global_age_pattern.csv")

# ==== 3. Attributable rate = attributable fraction x all cause mortality (every country) =========
cat("\n==== 3. Country level: age standardised attributable rate = attributable fraction x all cause rate ====\n")
W <- ends(S %>% mutate(hp = h_rate / m_rate, cp = c_rate / m_rate), c("h_rate", "c_rate", "m_rate", "hp", "cp", "h_n", "c_n", "d_n", "h_rate_lo", "h_rate_hi", "c_rate_lo", "c_rate_hi"))
stopifnot(nrow(W) == 204)
cat(sprintf("Non-positive point estimates at the period ends: heat rate 1990 %d, 2023 %d; cold rate 1990 %d, 2023 %d (countries: %s)\n", sum(W$h_rate_1990 <= 0), sum(W$h_rate_2023 <= 0), sum(W$c_rate_1990 <= 0), sum(W$c_rate_2023 <= 0),
            paste(W$location_name[W$h_rate_1990 <= 0 | W$h_rate_2023 <= 0 | W$c_rate_1990 <= 0 | W$c_rate_2023 <= 0], collapse = "; ")))
P <- W %>% filter(h_rate_1990 > 0, h_rate_2023 > 0, c_rate_1990 > 0, c_rate_2023 > 0) %>%
  mutate(dM = log(m_rate_2023 / m_rate_1990), dHp = log(hp_2023 / hp_1990), dCp = log(cp_2023 / cp_1990), dH = log(h_rate_2023 / h_rate_1990), dC = log(c_rate_2023 / c_rate_1990),
         z90 = log(h_rate_1990 / c_rate_1990), z23 = log(h_rate_2023 / c_rate_2023))
cat("Countries with all four rates positive:", nrow(P), "\n")
cat(sprintf("All cause age standardised rate: fell in %d of %d countries; median log change %.3f (%.1f%%)\n", sum(P$dM < 0), nrow(P), median(P$dM), 100 * (exp(median(P$dM)) - 1)))
cat(sprintf("Heat: attributable fraction rose in %d; rate rose in %d; median log change fraction %.3f (x%.2f), rate %.3f (x%.2f)\n", sum(P$dHp > 0), sum(P$dH > 0), median(P$dHp), exp(median(P$dHp)), median(P$dH), exp(median(P$dH))))
cat(sprintf("Heat: countries where the fraction rose but the rate fell: %d; countries where the fraction rose and the rate rose by less: %d\n", sum(P$dHp > 0 & P$dH < 0), sum(P$dHp > 0 & P$dH > 0 & P$dH < P$dHp)))
cat(sprintf("Cold: attributable fraction fell in %d; rate fell in %d; median log change fraction %.3f (%.1f%%), rate %.3f (%.1f%%)\n", sum(P$dCp < 0), sum(P$dC < 0), median(P$dCp), 100 * (exp(median(P$dCp)) - 1), median(P$dC), 100 * (exp(median(P$dC)) - 1)))
fc <- P %>% filter(dC < 0)
cat(sprintf("Cold: share of the decline in the rate carried by falling all cause mortality, median over %d countries with a decline: %.0f%% (IQR %.0f to %.0f); countries where mortality carries more than half: %d\n",
            nrow(fc), 100 * median(fc$dM / fc$dC), 100 * quantile(fc$dM / fc$dC, .25), 100 * quantile(fc$dM / fc$dC, .75), sum(fc$dM / fc$dC > 0.5)))
cat(sprintf("Global (age standardised): heat rate x%.3f = fraction x%.3f times all cause x%.3f; cold rate x%.3f = fraction x%.3f times all cause x%.3f\n",
            g23$h_asr / g90$h_asr, (g23$h_asr / g23$m_asr) / (g90$h_asr / g90$m_asr), g23$m_asr / g90$m_asr, g23$c_asr / g90$c_asr, (g23$c_asr / g23$m_asr) / (g90$c_asr / g90$m_asr), g23$m_asr / g90$m_asr))
cat(sprintf("Global: heat rate in 2023 had all cause mortality stayed at its 1990 level: %.2f per 100,000 (observed %.2f; 1990 %.2f); share of the fall in the cold rate carried by all cause mortality: %.0f%%\n",
            g90$h_asr * (g23$h_asr / g23$m_asr) / (g90$h_asr / g90$m_asr), g23$h_asr, g90$h_asr, 100 * log(g23$m_asr / g90$m_asr) / log(g23$c_asr / g90$c_asr)))
cat("By super region (median log changes; number of countries; countries with a rising heat rate):\n")
print(as.data.frame(P %>% group_by(region = substr(sr_short[super_region_name], 1, 30)) %>% summarise(n = n(), heat_rate_up = sum(dH > 0), heat_frac_up = sum(dHp > 0), cold_frac_down = sum(dCp < 0), dM = round(median(dM), 2), dHp = round(median(dHp), 2), dH = round(median(dH), 2), dCp = round(median(dCp), 2), dC = round(median(dC), 2))), row.names = FALSE)
write_csv(P, "output/country_identity.csv")

# ==== 4. Composition across countries ============================================================
cat("\n==== 4. Heat share and log ratio across countries ====\n")
C <- S %>% mutate(share = 100 * h_rate / (h_rate + c_rate), share_n = 100 * h_n / (h_n + c_n), pos = h_rate > 0 & c_rate > 0, z = ifelse(pos, log(pmax(h_rate, 1e-12) / pmax(c_rate, 1e-12)), NA_real_))
Cn <- C %>% filter(!shock)
cy <- Cn %>% group_by(year) %>% summarise(n = n(), heat_gt_cold = sum(h_rate > c_rate), med_share = median(share), q25 = quantile(share, .25), q75 = quantile(share, .75), .groups = "drop")
Wc <- W %>% mutate(share_1990 = 100 * h_rate_1990 / (h_rate_1990 + c_rate_1990), share_2023 = 100 * h_rate_2023 / (h_rate_2023 + c_rate_2023), dom90 = h_rate_1990 > c_rate_1990, dom23 = h_rate_2023 > c_rate_2023)
cat(sprintf("Countries where heat exceeds cold (age standardised rate): %d in 1990, %d in 2023 (of 204, period ends); by year excluding shock years: 2000 %d of %d, 2010 %d of %d, 2019 %d of %d\n", sum(Wc$dom90), sum(Wc$dom23),
            cy$heat_gt_cold[cy$year == 2000], cy$n[cy$year == 2000], cy$heat_gt_cold[cy$year == 2010], cy$n[cy$year == 2010], cy$heat_gt_cold[cy$year == 2019], cy$n[cy$year == 2019]))
cat(sprintf("Median country heat share: %.1f%% (IQR %.1f to %.1f) in 1990; %.1f%% (IQR %.1f to %.1f) in 2023\n", median(Wc$share_1990), quantile(Wc$share_1990, .25), quantile(Wc$share_1990, .75), median(Wc$share_2023), quantile(Wc$share_2023, .25), quantile(Wc$share_2023, .75)))
both <- Wc %>% filter(h_rate_1990 + c_rate_1990 > 0, h_rate_2023 + c_rate_2023 > 0)
cat(sprintf("Heat share rose between 1990 and 2023 in %d of %d countries with a positive total at both ends; crossed from cold dominated to heat dominated: %d; the reverse: %d\n", sum(both$share_2023 > both$share_1990), nrow(both), sum(Wc$dom23 & !Wc$dom90), sum(!Wc$dom23 & Wc$dom90)))
cat("Countries that crossed to heat dominated:", paste(sort(Wc$location_name[Wc$dom23 & !Wc$dom90]), collapse = "; "), "\n")
cat("Countries that crossed to cold dominated:", paste(sort(Wc$location_name[!Wc$dom23 & Wc$dom90]), collapse = "; "), "\n")
cat(sprintf("Share of heat deaths in 2023 occurring in heat dominated countries: %.1f%%; share of all deaths (2023) in heat dominated countries: %.1f%% (1990 classification: %.1f%%)\n", 100 * sum(Wc$h_n_2023[Wc$dom23]) / sum(Wc$h_n_2023),
            100 * sum(Wc$d_n_2023[Wc$dom23]) / sum(Wc$d_n_2023), 100 * sum(Wc$d_n_1990[Wc$dom90]) / sum(Wc$d_n_1990)))
# balanced set for dispersion: positive heat and cold rates in every year that is not a shock year; shock years interpolated
bal <- Cn %>% group_by(location_id) %>% summarise(ok = all(pos)) %>% filter(ok) %>% pull(location_id)
cat("Balanced set (positive heat and cold rates in all years other than shock years):", length(bal), "countries; excluded:", paste(ctry$location_name[!ctry$location_id %in% bal], collapse = "; "), "\n")
Z <- C %>% filter(location_id %in% bal) %>% arrange(location_id, year) %>% group_by(location_id) %>% mutate(zi = approx(year[!shock], z[!shock], xout = year, rule = 2)$y) %>% ungroup()
cat("Shock years interpolated in the balanced set:", sum(Z$shock), "\n")
sig <- Z %>% group_by(year) %>% summarise(sd_z = sd(zi), mean_z = mean(zi), med_z = median(zi), iqr_z = IQR(zi), .groups = "drop")
cat(sprintf("SD of log(heat/cold) across %d countries: 1990 %.3f, 2023 %.3f (change %+.1f%%); min %.3f (%d), max %.3f (%d); IQR 1990 %.3f, 2023 %.3f\n", length(bal), sig$sd_z[1], sig$sd_z[34], pct(sig$sd_z[1], sig$sd_z[34]), min(sig$sd_z), sig$year[which.min(sig$sd_z)], max(sig$sd_z), sig$year[which.max(sig$sd_z)], sig$iqr_z[1], sig$iqr_z[34]))
st <- summary(lm(sd_z ~ year, sig))$coefficients; cat(sprintf("Linear trend in the SD: %.4f per year (SE %.4f, p = %.3f)\n", st[2, 1], st[2, 2], st[2, 4]))
B <- Z %>% filter(year %in% c(1990, 2023)) %>% select(location_id, super_region_name, year, zi) %>% pivot_wider(names_from = year, values_from = zi, names_prefix = "z")
cat(sprintf("Log ratio: mean 1990 %.3f, 2023 %.3f; median change %+.3f (heat to cold ratio multiplied by %.2f in the median country), IQR of change %+.3f to %+.3f; change positive in %d of %d\n", mean(B$z1990), mean(B$z2023), median(B$z2023 - B$z1990), exp(median(B$z2023 - B$z1990)),
            quantile(B$z2023 - B$z1990, .25), quantile(B$z2023 - B$z1990, .75), sum(B$z2023 > B$z1990), nrow(B)))
bfit <- summary(lm(I(z2023 - z1990) ~ z1990, B))
cat(sprintf("Beta: change in log ratio on 1990 log ratio: slope %.4f (SE %.4f), p = %.3f, r = %.3f, n = %d\n", bfit$coefficients[2, 1], bfit$coefficients[2, 2], bfit$coefficients[2, 4], cor(B$z2023 - B$z1990, B$z1990), nrow(B)))
an <- function(v) { a <- anova(lm(v ~ super_region_name, B)); 100 * a$`Sum Sq`[1] / sum(a$`Sum Sq`) }
cat(sprintf("Share of the variance of the log ratio between super regions: %.1f%% (1990), %.1f%% (2023); rank correlation of country log ratios 1990 vs 2023: %.3f\n", an(B$z1990), an(B$z2023), cor(B$z1990, B$z2023, method = "spearman")))
write_csv(sig, "output/sigma_series.csv"); write_csv(cy, "output/share_by_year.csv")
write_csv(C %>% select(location_id, location_name, super_region_name, year, shock, h_rate, c_rate, share, share_n, z, h_n, c_n, h_paf, c_paf, m_rate, d_n), "output/country_year.csv")

# climate axis: log ratio against population weighted mean temperature (covariate valid to 2020)
per <- function(y) case_when(y <= 1994 ~ "p1", y >= 2016 & y <= 2020 ~ "p2", TRUE ~ NA_character_)
T5 <- temp %>% mutate(period = per(year)) %>% filter(!is.na(period)) %>% group_by(location_id, period) %>% summarise(temp = mean(temp), .groups = "drop")
Z5 <- Cn %>% filter(location_id %in% bal) %>% mutate(period = per(year)) %>% filter(!is.na(period)) %>% group_by(location_id, period) %>% summarise(z = mean(z), .groups = "drop")
ZT <- inner_join(Z5, T5, by = c("location_id", "period")) %>% left_join(ctry, by = "location_id")
m1 <- lm(z ~ temp, ZT %>% filter(period == "p1")); m2 <- lm(z ~ temp, ZT %>% filter(period == "p2")); mp <- lm(z ~ temp + period, ZT)
wide <- ZT %>% select(location_id, period, z, temp) %>% pivot_wider(names_from = period, values_from = c(z, temp))
cat(sprintf("\nClimate axis (five year means, 1990-94 and 2016-20, n = %d): log ratio on mean temperature, slope %.3f per degree (R2 %.3f) in 1990-94; %.3f (R2 %.3f) in 2016-20\n", nrow(wide), coef(m1)[2], summary(m1)$r.squared, coef(m2)[2], summary(m2)$r.squared))
cat(sprintf("Mean temperature change between periods: %.3f degrees (IQR %.3f to %.3f); warmed in %d of %d countries\n", mean(wide$temp_p2 - wide$temp_p1), quantile(wide$temp_p2 - wide$temp_p1, .25), quantile(wide$temp_p2 - wide$temp_p1, .75), sum(wide$temp_p2 > wide$temp_p1), nrow(wide)))
cat(sprintf("Pooled model with common slope %.3f: mean change in log ratio %.3f = movement along the cross sectional gradient %.3f (%.0f%%) + upward shift of the gradient %.3f (%.0f%%)\n", coef(mp)["temp"], mean(wide$z_p2 - wide$z_p1), coef(mp)["temp"] * mean(wide$temp_p2 - wide$temp_p1),
            100 * coef(mp)["temp"] * mean(wide$temp_p2 - wide$temp_p1) / mean(wide$z_p2 - wide$z_p1), coef(mp)["periodp2"], 100 * coef(mp)["periodp2"] / mean(wide$z_p2 - wide$z_p1)))
write_csv(ZT, "output/climate_axis.csv")

# within country response of the attributable fractions to annual temperature, 1990 to 2020
FE <- Cn %>% filter(year <= 2020, location_id %in% bal) %>% mutate(lhp = log(h_rate / m_rate), lcp = log(c_rate / m_rate), z = lhp - lcp) %>% inner_join(temp %>% select(location_id, year, temp), by = c("location_id", "year"))
# residualise on country intercepts and country specific linear trends (Frisch-Waugh), then pool
res <- FE %>% group_by(location_id) %>% mutate(across(c(lhp, lcp, z, temp), ~ resid(lm(.x ~ year)), .names = "r_{.col}")) %>% ungroup()
slope <- function(d, y) sum(d$r_temp * d[[y]]) / sum(d$r_temp^2)
ids <- unique(res$location_id); spl <- split(res, res$location_id)
set.seed(11); boot <- t(replicate(2000, { d <- bind_rows(spl[sample(length(spl), replace = TRUE)]); c(slope(d, "r_lhp"), slope(d, "r_lcp"), slope(d, "r_z")) }))
bh <- slope(res, "r_lhp"); bc <- slope(res, "r_lcp"); bz <- slope(res, "r_z")
cat(sprintf("\nWithin country annual response (country intercepts and country specific linear trends, 1990 to 2020, %d countries, %d observations; intervals from 2000 bootstrap resamples of countries):\n", length(ids), nrow(res)))
cat(sprintf("  per degree Celsius warmer year: heat fraction %+.1f%% (95%% CI %+.1f to %+.1f); cold fraction %+.1f%% (%+.1f to %+.1f); log ratio %+.3f (%+.3f to %+.3f)\n",
            100 * (exp(bh) - 1), 100 * (exp(quantile(boot[, 1], .025)) - 1), 100 * (exp(quantile(boot[, 1], .975)) - 1), 100 * (exp(bc) - 1), 100 * (exp(quantile(boot[, 2], .025)) - 1), 100 * (exp(quantile(boot[, 2], .975)) - 1), bz, quantile(boot[, 3], .025), quantile(boot[, 3], .975)))
cs <- res %>% group_by(location_id) %>% summarise(b_h = sum(r_temp * r_lhp) / sum(r_temp^2), b_c = sum(r_temp * r_lcp) / sum(r_temp^2), b_z = sum(r_temp * r_z) / sum(r_temp^2), sd_temp = sd(r_temp), .groups = "drop")
cat(sprintf("  country specific slopes: heat fraction median %+.3f log points per degree (IQR %+.3f to %+.3f), positive in %d of %d; cold fraction median %+.3f (IQR %+.3f to %+.3f), negative in %d\n", median(cs$b_h), quantile(cs$b_h, .25), quantile(cs$b_h, .75), sum(cs$b_h > 0), nrow(cs),
            median(cs$b_c), quantile(cs$b_c, .25), quantile(cs$b_c, .75), sum(cs$b_c < 0)))
tr <- FE %>% group_by(location_id) %>% summarise(dT = coef(lm(temp ~ year))[2] * 30, dz = coef(lm(z ~ year))[2] * 30, dlhp = coef(lm(lhp ~ year))[2] * 30, dlcp = coef(lm(lcp ~ year))[2] * 30, .groups = "drop") %>% left_join(cs, by = "location_id")
cat(sprintf("  linear trends over 1990 to 2020 (30 years): temperature median %+.2f degrees (IQR %+.2f to %+.2f), positive in %d of %d countries; log heat fraction median %+.3f; log cold fraction median %+.3f; log ratio median %+.3f\n",
            median(tr$dT), quantile(tr$dT, .25), quantile(tr$dT, .75), sum(tr$dT > 0), nrow(tr), median(tr$dlhp), median(tr$dlcp), median(tr$dz)))
cat(sprintf("  implied by the pooled annual response at the median warming: heat fraction %+.3f (%.0f%% of the median trend), cold fraction %+.3f (%.0f%%), log ratio %+.3f (%.0f%%)\n",
            bh * median(tr$dT), 100 * bh * median(tr$dT) / median(tr$dlhp), bc * median(tr$dT), 100 * bc * median(tr$dT) / median(tr$dlcp), bz * median(tr$dT), 100 * bz * median(tr$dT) / median(tr$dz)))
tr <- tr %>% mutate(imp_z = b_z * dT, imp_h = b_h * dT, imp_c = b_c * dT)
cat(sprintf("  country by country (own slope x own warming): median implied change in log ratio %+.3f against median observed %+.3f; median implied share of the observed trend %.0f%% (IQR %.0f to %.0f) among %d countries with a positive trend in the log ratio\n",
            median(tr$imp_z), median(tr$dz), 100 * median((tr$imp_z / tr$dz)[tr$dz > 0]), 100 * quantile((tr$imp_z / tr$dz)[tr$dz > 0], .25), 100 * quantile((tr$imp_z / tr$dz)[tr$dz > 0], .75), sum(tr$dz > 0)))
cat(sprintf("  across countries: Spearman correlation between the temperature trend and the trend in the log ratio %.3f; between implied and observed change in the log ratio %.3f; in log heat fraction %.3f; in log cold fraction %.3f\n",
            cor(tr$dT, tr$dz, method = "spearman"), cor(tr$imp_z, tr$dz, method = "spearman"), cor(tr$imp_h, tr$dlhp, method = "spearman"), cor(tr$imp_c, tr$dlcp, method = "spearman")))
write_csv(tr %>% left_join(ctry, by = "location_id"), "output/country_trends.csv")

# ==== 5. Cause composition ========================================================================
cat("\n==== 5. Causes: within cause change against change in the cause mix (all ages, 1990 and 2023) ====\n")
tot <- S %>% filter(year %in% c(1990, 2023)) %>% select(location_id, year, H = h_n, Cc = c_n, D = d_n)
K <- causes %>% inner_join(tot, by = c("location_id", "year"))
ck <- K %>% group_by(location_id, year) %>% summarise(dh = sum(h_n) - first(H), dc = sum(c_n) - first(Cc), H = first(H), Cc = first(Cc), .groups = "drop")
cat(sprintf("Sum of the seven cause specific attributable deaths against the all cause attributable total: largest absolute difference heat %.2f deaths, cold %.2f deaths\n", max(abs(ck$dh)), max(abs(ck$dc))))
gk <- K %>% group_by(year, cause) %>% summarise(h = sum(h_n), c = sum(c_n), d = sum(d_n), .groups = "drop") %>% left_join(tot %>% group_by(year) %>% summarise(D = sum(D), .groups = "drop"), by = "year") %>% mutate(w = d / D, ph = h / d, pc = c / d)
gw <- gk %>% pivot_wider(names_from = year, values_from = c(h, c, d, D, w, ph, pc))
kit <- function(p0, p1, w0, w1) tibble(within = (p1 - p0) * (w0 + w1) / 2, mix = (w1 - w0) * (p0 + p1) / 2)
gh <- bind_cols(gw %>% select(cause), kit(gw$ph_1990, gw$ph_2023, gw$w_1990, gw$w_2023)); gc <- bind_cols(gw %>% select(cause), kit(gw$pc_1990, gw$pc_2023, gw$w_1990, gw$w_2023))
cat("Global (sum of 204 countries): cause share of all deaths w (%), heat and cold deaths, cause specific attributable fractions (%):\n")
print(as.data.frame(gw %>% transmute(cause = substr(cause, 1, 40), w90 = round(100 * w_1990, 1), w23 = round(100 * w_2023, 1), heat90 = round(h_1990), heat23 = round(h_2023), cold90 = round(c_1990), cold23 = round(c_2023),
                       ph90 = round(100 * ph_1990, 2), ph23 = round(100 * ph_2023, 2), pc90 = round(100 * pc_1990, 2), pc23 = round(100 * pc_2023, 2))), row.names = FALSE)
cat(sprintf("Seven outcome groups as a share of all deaths: %.1f%% (1990), %.1f%% (2023)\n", 100 * sum(gw$w_1990), 100 * sum(gw$w_2023)))
cat(sprintf("Gross and net cold: sum of positive cause components %s (1990), %s (2023); sum of negative (injury) components %s (1990), %s (2023)\n", fmt(sum(gw$c_1990[gw$c_1990 > 0]), 0), fmt(sum(gw$c_2023[gw$c_2023 > 0]), 0), fmt(sum(gw$c_1990[gw$c_1990 < 0]), 0), fmt(sum(gw$c_2023[gw$c_2023 < 0]), 0)))
cat(sprintf("Global heat fraction (crude): %.3f%% to %.3f%%, change %.3f points = within cause %.3f (%.0f%%) + cause mix %.3f (%.0f%%)\n", 100 * sum(gw$w_1990 * gw$ph_1990), 100 * sum(gw$w_2023 * gw$ph_2023),
            100 * sum(gh$within + gh$mix), 100 * sum(gh$within), 100 * sum(gh$within) / sum(gh$within + gh$mix), 100 * sum(gh$mix), 100 * sum(gh$mix) / sum(gh$within + gh$mix)))
cat(sprintf("Global cold fraction (crude): %.3f%% to %.3f%%, change %+.3f points = within cause %+.3f + cause mix %+.3f\n", 100 * sum(gw$w_1990 * gw$pc_1990), 100 * sum(gw$w_2023 * gw$pc_2023),
            100 * sum(gc$within + gc$mix), 100 * sum(gc$within), 100 * sum(gc$mix)))
cat("Contribution by cause to the change in the global heat and cold fractions (hundredths of a percentage point):\n")
print(as.data.frame(tibble(cause = substr(gh$cause, 1, 40), heat_within = round(1e4 * gh$within, 2), heat_mix = round(1e4 * gh$mix, 2), cold_within = round(1e4 * gc$within, 2), cold_mix = round(1e4 * gc$mix, 2))), row.names = FALSE)
flag_end <- ep$location_id[ep$y0 != 1990 | ep$y1 != 2023]
kc <- K %>% filter(!location_id %in% flag_end) %>% mutate(w = d_n / D, ph = h_n / d_n, pc = c_n / d_n) %>% select(location_id, year, cause, w, ph, pc) %>% pivot_wider(names_from = year, values_from = c(w, ph, pc)) %>%
  group_by(location_id) %>% summarise(h_within = sum((ph_2023 - ph_1990) * (w_1990 + w_2023) / 2), h_mix = sum((w_2023 - w_1990) * (ph_1990 + ph_2023) / 2),
                                      c_within = sum((pc_2023 - pc_1990) * (w_1990 + w_2023) / 2), c_mix = sum((w_2023 - w_1990) * (pc_1990 + pc_2023) / 2), .groups = "drop")
cat(sprintf("Countries (n = %d, excluding those with a shock year at either end): heat fraction (crude) rose in %d; within cause component positive in %d; cause mix component positive in %d; mix larger than within in absolute value in %d\n", nrow(kc), sum(kc$h_within + kc$h_mix > 0), sum(kc$h_within > 0), sum(kc$h_mix > 0), sum(abs(kc$h_mix) > abs(kc$h_within))))
cat(sprintf("Countries: cold fraction (crude) fell in %d; within cause component negative in %d; cause mix component positive in %d; countries where the within cause component is negative but the crude fraction rose: %d\n", sum(kc$c_within + kc$c_mix < 0), sum(kc$c_within < 0), sum(kc$c_mix > 0), sum(kc$c_within < 0 & kc$c_within + kc$c_mix > 0)))
write_csv(kc %>% left_join(ctry, by = "location_id"), "output/country_cause_kitagawa.csv")
# Shapley attribution of the change in the global heat share to 14 components (7 causes x heat, cold), crude rates
N90 <- sum(byage$pop[byage$location_id == 1 & byage$year == 1990]); N23 <- sum(byage$pop[byage$location_id == 1 & byage$year == 2023])
x0 <- as.list(c(gw$h_1990, gw$c_1990) / N90); x1 <- as.list(c(gw$h_2023, gw$c_2023) / N23); names(x0) <- names(x1) <- c(paste0("heat: ", gw$cause), paste0("cold: ", gw$cause))
fS <- function(x) { v <- unlist(x); 100 * sum(v[1:7]) / sum(v) }
ph <- shapley(fS, x0, x1)
cat(sprintf("\nGlobal heat share (crude rates, sum of countries): %.2f%% to %.2f%%; Shapley attribution of the %.2f point change over 14 components (population %s to %s):\n", fS(x0), fS(x1), fS(x1) - fS(x0), fmt(N90, 0), fmt(N23, 0)))
sh <- tibble(component = names(ph), points = ph) %>% mutate(side = sub(":.*", "", component), cause = sub(".*: ", "", component)); print(as.data.frame(sh %>% mutate(points = round(points, 2), cause = substr(cause, 1, 40)) %>% select(side, cause, points)), row.names = FALSE)
cat(sprintf("  heat components total %+.2f points; cold components total %+.2f points\n", sum(sh$points[sh$side == "heat"]), sum(sh$points[sh$side == "cold"])))
write_csv(sh, "output/shapley_share_causes.csv"); write_csv(gw, "output/global_causes.csv")
gs <- gk %>% group_by(year) %>% mutate(heat_cause_share = 100 * h / sum(h), cold_cause_share = 100 * c / sum(c[c > 0]), heat_share_in_cause = ifelse(c > 0, 100 * h / (h + c), NA)) %>% ungroup()
cat("\nCause composition of heat deaths (% of heat deaths), of gross cold deaths (% of the positive components) and heat share within cause (%):\n"); print(as.data.frame(gs %>% transmute(cause = substr(cause, 1, 40), year, heat_cause_share = round(heat_cause_share, 1), cold_cause_share = round(cold_cause_share, 1), heat_share_in_cause = round(heat_share_in_cause, 1)) %>% arrange(cause, year)), row.names = FALSE)

# ==== 6. Uncertainty ==============================================================================
cat("\n==== 6. Uncertainty of the components and of derived indicators ====\n")
PA <- parent %>% filter(sex == "combined") %>% select(location_id, year, p_rate_pt, p_rate_lo, p_rate_hi)
u_end <- function(which) { yv <- if (which == 0) ep$y0 else ep$y1; nom <- if (which == 0) 1990 else 2023
  S %>% inner_join(tibble(location_id = ep$location_id, yy = yv), by = "location_id") %>% filter(year == yy) %>% left_join(PA %>% filter(year == nom) %>% select(-year), by = "location_id") %>%
    mutate(sub = yy != nom, sh = (h_rate_hi - h_rate_lo) / 3.92, sc = (c_rate_hi - c_rate_lo) / 3.92, sp = (p_rate_hi - p_rate_lo) / 3.92, rho = ifelse(sub, NA, (sp^2 - sh^2 - sc^2) / (2 * sh * sc)),
           rw_h = (h_rate_hi - h_rate_lo) / abs(h_rate), rw_c = (c_rate_hi - c_rate_lo) / abs(c_rate), add = (h_rate + c_rate) / p_rate_pt) }
u90 <- u_end(0); u23 <- u_end(1)
cat(sprintf("Additivity of heat and cold to the parent risk (age standardised rate, 2023): median %.4f, 1st to 99th percentile %.4f to %.4f\n", median(u23$add, na.rm = TRUE), quantile(u23$add, .01, na.rm = TRUE), quantile(u23$add, .99, na.rm = TRUE)))
cat(sprintf("Relative width of the 95%% UI of the age standardised rate, 2023: heat median %.2f (IQR %.2f to %.2f); cold median %.2f (IQR %.2f to %.2f)\n", median(u23$rw_h), quantile(u23$rw_h, .25), quantile(u23$rw_h, .75), median(u23$rw_c), quantile(u23$rw_c, .25), quantile(u23$rw_c, .75)))
cat(sprintf("Countries whose heat interval includes zero: %d (1990), %d (2023); cold interval includes zero: %d (1990), %d (2023)\n", sum(u90$h_rate_lo <= 0), sum(u23$h_rate_lo <= 0), sum(u90$c_rate_lo <= 0), sum(u23$c_rate_lo <= 0)))
cat(sprintf("Implied correlation between heat and cold errors from the parent interval, 2023: median %.2f (IQR %.2f to %.2f); share negative %.0f%%; share outside [-1, 1] %.0f%%; 1990 median %.2f\n", median(u23$rho, na.rm = TRUE), quantile(u23$rho, .25, na.rm = TRUE), quantile(u23$rho, .75, na.rm = TRUE),
            100 * mean(u23$rho < 0, na.rm = TRUE), 100 * mean(abs(u23$rho) > 1, na.rm = TRUE), median(u90$rho, na.rm = TRUE)))
u90$rho[!is.finite(u90$rho)] <- median(u90$rho, na.rm = TRUE); u23$rho[!is.finite(u23$rho)] <- median(u23$rho, na.rm = TRUE)
gpp <- agg %>% filter(location_id == 1, year %in% c(1990, 2023)) %>% mutate(sh = (h_n_hi - h_n_lo) / 3.92, sc = (c_n_hi - c_n_lo) / 3.92, sp = (p_n_hi - p_n_lo) / 3.92, rho = (sp^2 - sh^2 - sc^2) / (2 * sh * sc))
cat(sprintf("Implied correlation at the global level (counts): %.2f (1990), %.2f (2023)\n", gpp$rho[1], gpp$rho[2]))
clip <- function(r) max(-0.9, min(0.9, r))
mc_share <- function(z, h, hl, hu, c, cl, cu) { H <- tpn(z[, 1], h, hl, hu); K <- tpn(z[, 2], c, cl, cu); 100 * H / (H + K) }
set.seed(1); e1 <- rnorm(NDRAW); e2 <- rnorm(NDRAW); mkz <- function(r) cbind(e1, clip(r) * e1 + sqrt(1 - clip(r)^2) * e2)
s90 <- mc_share(mkz(gpp$rho[1]), g90$h_n, g90$h_n_lo, g90$h_n_hi, g90$c_n, g90$c_n_lo, g90$c_n_hi)
s23 <- mc_share(mkz(gpp$rho[2]), g23$h_n, g23$h_n_lo, g23$h_n_hi, g23$c_n, g23$c_n_lo, g23$c_n_hi)
s23i <- mc_share(rcop(NDRAW, clip(gpp$rho[2])), g23$h_n, g23$h_n_lo, g23$h_n_hi, g23$c_n, g23$c_n_lo, g23$c_n_hi)
cat(sprintf("Global heat share, approximate 95%% interval: 1990 %.1f%% (%.1f to %.1f); 2023 %.1f%% (%.1f to %.1f); change %.1f points (%.1f to %.1f) with errors shared across years; (%.1f to %.1f) with independent years\n",
            g90$share_n, qi(s90)[1], qi(s90)[2], g23$share_n, qi(s23)[1], qi(s23)[2], g23$share_n - g90$share_n, qi(s23 - s90)[1], qi(s23 - s90)[2], qi(s23i - s90)[1], qi(s23i - s90)[2]))
# country shares at the two ends, with intervals; deviates shared across years, or independent
set.seed(3)
cmc <- map_dfr(u23$location_id, function(L) { a <- u90[u90$location_id == L, ]; b <- u23[u23$location_id == L, ]
  e1 <- rnorm(NDRAW); e2 <- rnorm(NDRAW); mk <- function(r) cbind(e1, clip(r) * e1 + sqrt(1 - clip(r)^2) * e2); z1 <- mk(b$rho); z0 <- mk(a$rho); z0i <- rcop(NDRAW, clip(a$rho))
  H1 <- tpn(z1[, 1], b$h_rate, b$h_rate_lo, b$h_rate_hi); K1 <- tpn(z1[, 2], b$c_rate, b$c_rate_lo, b$c_rate_hi)
  H0 <- tpn(z0[, 1], a$h_rate, a$h_rate_lo, a$h_rate_hi); K0 <- tpn(z0[, 2], a$c_rate, a$c_rate_lo, a$c_rate_hi)
  H0i <- tpn(z0i[, 1], a$h_rate, a$h_rate_lo, a$h_rate_hi); K0i <- tpn(z0i[, 2], a$c_rate, a$c_rate_lo, a$c_rate_hi)
  s1 <- 100 * H1 / (H1 + K1); s0 <- 100 * H0 / (H0 + K0)
  tibble(location_id = L, share90 = 100 * a$h_rate / (a$h_rate + a$c_rate), share23 = 100 * b$h_rate / (b$h_rate + b$c_rate), lo23 = qi(s1)[1], hi23 = qi(s1)[2], lo90 = qi(s0)[1], hi90 = qi(s0)[2],
         p_rise_shared = mean(H1 * K0 > H0 * K1), p_rise_indep = mean(H1 * K0i > H0i * K1), p_heat_gt_cold23 = mean(H1 > K1), p_heat_gt_cold90 = mean(H0 > K0)) })
cat(sprintf("Countries where the heat to cold ratio in 2023 exceeds that of 1990 in at least 97.5%% of draws: %d (errors shared across years), %d (independent years); in at most 2.5%%: %d, %d\n", sum(cmc$p_rise_shared >= 0.975), sum(cmc$p_rise_indep >= 0.975), sum(cmc$p_rise_shared <= 0.025), sum(cmc$p_rise_indep <= 0.025)))
cat(sprintf("Countries classified with at least 97.5%% probability: 2023 heat dominated %d, cold dominated %d, unresolved %d; 1990 heat dominated %d, cold dominated %d, unresolved %d\n", sum(cmc$p_heat_gt_cold23 >= 0.975), sum(cmc$p_heat_gt_cold23 <= 0.025), sum(cmc$p_heat_gt_cold23 > 0.025 & cmc$p_heat_gt_cold23 < 0.975),
            sum(cmc$p_heat_gt_cold90 >= 0.975), sum(cmc$p_heat_gt_cold90 <= 0.025), sum(cmc$p_heat_gt_cold90 > 0.025 & cmc$p_heat_gt_cold90 < 0.975)))
write_csv(cmc %>% left_join(ctry, by = "location_id"), "output/country_share_mc.csv")
# dispersion and beta with uncertainty: countries with strictly positive lower bounds for both risks at both ends
idp <- intersect(u90$location_id[u90$h_rate_lo > 0 & u90$c_rate_lo > 0], u23$location_id[u23$h_rate_lo > 0 & u23$c_rate_lo > 0]); cat("Countries with positive lower bounds for heat and cold at both ends:", length(idp), "\n")
e0 <- u90[match(idp, u90$location_id), ]; e1 <- u23[match(idp, u23$location_id), ]
set.seed(4); z1l <- lapply(idp, function(L) rcop(NDRAW, 0)); set.seed(5); z2l <- lapply(idp, function(L) rcop(NDRAW, 0))
dr <- function(e, zl) sapply(seq_along(idp), function(i) log(sln(zl[[i]][, 1], e$h_rate[i], e$h_rate_lo[i], e$h_rate_hi[i])) - log(sln(zl[[i]][, 2], e$c_rate[i], e$c_rate_lo[i], e$c_rate_hi[i])))
Z0 <- dr(e0, z1l); Z1 <- dr(e1, z1l); Z0i <- dr(e0, z2l)
sd0 <- apply(Z0, 1, sd); sd1 <- apply(Z1, 1, sd); sd0i <- apply(Z0i, 1, sd)
bet <- sapply(1:NDRAW, function(k) coef(lm(I(Z1[k, ] - Z0[k, ]) ~ Z0[k, ]))[2]); beti <- sapply(1:NDRAW, function(k) coef(lm(I(Z1[k, ] - Z0i[k, ]) ~ Z0i[k, ]))[2])
zp0 <- log(e0$h_rate / e0$c_rate); zp1 <- log(e1$h_rate / e1$c_rate)
cat(sprintf("In this set (n = %d): SD of log ratio %.3f (1990), %.3f (2023); share of draws in which the 2023 SD is below the 1990 SD: %.0f%% (errors shared across years), %.0f%% (independent years); interval for the change in SD %.3f to %.3f (shared), %.3f to %.3f (independent)\n",
            length(idp), sd(zp0), sd(zp1), 100 * mean(sd1 < sd0), 100 * mean(sd1 < sd0i), qi(sd1 - sd0)[1], qi(sd1 - sd0)[2], qi(sd1 - sd0i)[1], qi(sd1 - sd0i)[2]))
cat(sprintf("Beta slope in this set: %.4f; Monte Carlo interval %.4f to %.4f (shared), %.4f to %.4f (independent years)\n", coef(lm(I(zp1 - zp0) ~ zp0))[2], qi(bet)[1], qi(bet)[2], qi(beti)[1], qi(beti)[2]))

# ==== 7. Anomalous heat years ======================================================================
cat("\n==== 7. Years of anomalous heat attributable fraction ====\n")
AN <- Cn %>% mutate(hp = h_rate / m_rate) %>% group_by(location_id) %>% filter(all(hp > 0)) %>% ungroup() %>% mutate(lhp = log(hp), lhr = log(h_rate), fit_yr = year <= 2019)
AN <- AN %>% group_by(location_id) %>% group_modify(~ { f <- lm(lhp ~ year, data = .x[.x$fit_yr, ]); r <- .x$lhp - predict(f, .x); s <- 1.4826 * median(abs(r[.x$fit_yr] - median(r[.x$fit_yr])))
  f2 <- lm(lhr ~ year, data = .x[.x$fit_yr, ]); r2 <- .x$lhr - predict(f2, .x)
  mutate(.x, res = r, zs = r / s, exp_paf = exp(predict(f, .x)), excess = (1 - exp_paf / hp) * h_n, res_rate = r2) }) %>% ungroup()
cat("Countries with a positive heat fraction in all years other than shock years:", n_distinct(AN$location_id), "; log linear trend fitted to 1990 to 2019; anomaly = residual / robust SD (1.4826 x MAD)\n")
ay <- AN %>% group_by(year) %>% summarise(n = n(), n_hot = sum(zs > 2), n_cool = sum(zs < -2), med_res = median(res), med_res_rate = median(res_rate), excess = sum(excess), heat = sum(h_n), .groups = "drop")
cat("Years 1990 to 2019 with the most countries in positive anomaly (z > 2), with the sum over countries of heat deaths above trend:\n"); print(as.data.frame(ay %>% filter(year <= 2019) %>% arrange(-n_hot) %>% head(8) %>% mutate(med_res = round(med_res, 3), med_res_rate = round(med_res_rate, 3), excess = round(excess), heat = round(heat))), row.names = FALSE)
cat("2020 to 2023 against the 1990 to 2019 trend (median residual of the log heat fraction; of the log heat rate):\n"); print(as.data.frame(ay %>% filter(year >= 2020) %>% transmute(year, n, med_res = round(med_res, 3), med_res_rate = round(med_res_rate, 3))), row.names = FALSE)
for (y in 2020:2023) { a <- AN %>% filter(year == y); cat(sprintf("  %d: heat fraction above trend in %d of %d countries; heat rate above trend in %d\n", y, sum(a$res > 0), nrow(a), sum(a$res_rate > 0))) }
top <- AN %>% filter(zs > 2, year <= 2019) %>% arrange(-excess) %>% head(15) %>% transmute(location_name = substr(location_name, 1, 28), year, z = round(zs, 1), heat_deaths = round(h_n), above_trend = round(excess), ratio_to_trend = round(hp / exp_paf, 2))
cat("Largest country year anomalies (z > 2) by heat deaths above trend, 1990 to 2019:\n"); print(as.data.frame(top), row.names = FALSE)
for (cn in c("Russian Federation", "France", "Italy", "Spain", "Germany", "India", "Pakistan", "Australia", "United States of America", "China", "Japan", "Brazil")) { a <- AN %>% filter(location_name == cn) %>% arrange(-zs) %>% head(3)
  cat(sprintf("  %s: top anomaly years %s\n", cn, paste(sprintf("%d (z = %.1f, x%.2f)", a$year, a$zs, a$hp / a$exp_paf), collapse = ", "))) }
RT <- AN %>% filter(year <= 2019) %>% inner_join(temp %>% select(location_id, year, temp), by = c("location_id", "year")) %>% group_by(location_id) %>% summarise(r = cor(res, resid(lm(temp ~ year))), .groups = "drop")
cat(sprintf("Within country correlation between the residual of the log heat fraction and the detrended annual mean temperature, 1990 to 2019: median %.2f (IQR %.2f to %.2f); positive in %d of %d countries; above 0.5 in %d\n", median(RT$r), quantile(RT$r, .25), quantile(RT$r, .75), sum(RT$r > 0), nrow(RT), sum(RT$r > 0.5)))
write_csv(AN %>% select(location_id, location_name, super_region_name, year, hp, lhp, res, zs, exp_paf, excess, h_n, res_rate), "output/anomalies.csv"); write_csv(ay, "output/anomaly_years.csv")

# ==== 8. Sex =======================================================================================
cat("\n==== 8. Sex ====\n")
sxe <- function(sx) ends(series %>% filter(sex == sx) %>% mutate(hp = h_rate / m_rate, cp = c_rate / m_rate), c("h_rate", "c_rate", "h_n", "c_n", "hp", "cp")) %>% rename_with(~ paste0(.x, "_", sx), -c(location_id, location_name, super_region_name))
SX <- inner_join(sxe("male"), sxe("female") %>% select(-location_name, -super_region_name), by = "location_id")
for (y in c("1990", "2023")) { v <- function(nm, sx) SX[[paste0(nm, "_", y, "_", sx)]]
  sm <- 100 * v("h_rate", "male") / (v("h_rate", "male") + v("c_rate", "male")); sf <- 100 * v("h_rate", "female") / (v("h_rate", "female") + v("c_rate", "female")); okp <- v("hp", "female") > 0 & v("hp", "male") > 0; okc <- v("cp", "female") > 0 & v("cp", "male") > 0
  cat(sprintf("%s: global heat share of counts male %.1f%%, female %.1f%%; male heat share above female in %d of 204 countries; median male to female ratio of the age standardised rate: heat %.2f, cold %.2f; of the age standardised attributable fraction: heat %.2f, cold %.2f\n", y,
              100 * sum(v("h_n", "male")) / sum(v("h_n", "male") + v("c_n", "male")), 100 * sum(v("h_n", "female")) / sum(v("h_n", "female") + v("c_n", "female")), sum(sm > sf, na.rm = TRUE),
              median((v("h_rate", "male") / v("h_rate", "female"))[okp]), median((v("c_rate", "male") / v("c_rate", "female"))[okc]), median((v("hp", "male") / v("hp", "female"))[okp]), median((v("cp", "male") / v("cp", "female"))[okc]))) }
cat(sprintf("2023 global counts: heat deaths male %s, female %s; cold deaths male %s, female %s\n", fmt(sum(SX$h_n_2023_male), 0), fmt(sum(SX$h_n_2023_female), 0), fmt(sum(SX$c_n_2023_male), 0), fmt(sum(SX$c_n_2023_female), 0)))

# ==== 9. Years of life lost per death ================================================================
cat("\n==== 9. Years of life lost per attributable death (1990 and 2023, sum over countries) ====\n")
Y <- yll %>% inner_join(S %>% select(location_id, year, h_n, c_n, super_region_name), by = c("location_id", "year"))
yy <- Y %>% group_by(year) %>% summarise(h = sum(h_yll) / sum(h_n), c = sum(c_yll) / sum(c_n), h_yll = sum(h_yll), c_yll = sum(c_yll), .groups = "drop")
cat(sprintf("Global YLLs per death: heat %.1f (1990), %.1f (2023); cold %.1f (1990), %.1f (2023); heat share of temperature attributable YLLs %.1f%% (1990), %.1f%% (2023); YLLs heat %s to %s, cold %s to %s\n", yy$h[1], yy$h[2], yy$c[1], yy$c[2],
            100 * yy$h_yll[1] / (yy$h_yll[1] + yy$c_yll[1]), 100 * yy$h_yll[2] / (yy$h_yll[2] + yy$c_yll[2]), fmt(yy$h_yll[1], 0), fmt(yy$h_yll[2], 0), fmt(yy$c_yll[1], 0), fmt(yy$c_yll[2], 0)))
print(as.data.frame(Y %>% group_by(region = substr(sr_short[super_region_name], 1, 30), year) %>% summarise(heat = round(sum(h_yll) / sum(h_n), 1), cold = round(sum(c_yll) / sum(c_n), 1), heat_share_yll = round(100 * sum(h_yll) / sum(h_yll + c_yll), 1), .groups = "drop")), row.names = FALSE)

# ==== 10. Concentration ==============================================================================
cat("\n==== 10. Concentration across countries ====\n")
for (y in c("1990", "2023")) { hh <- pmax(W[[paste0("h_n_", y)]], 0); cc <- pmax(W[[paste0("c_n_", y)]], 0); dd <- W[[paste0("d_n_", y)]]
  cat(sprintf("%s: Gini of heat deaths %.3f, cold deaths %.3f, all cause deaths %.3f; five largest contributors to heat deaths %.1f%% (%s); to cold deaths %.1f%% (%s)\n", y, gini(hh), gini(cc), gini(dd),
              100 * sum(sort(hh, TRUE)[1:5]) / sum(hh), paste(W$location_name[order(-hh)][1:5], collapse = ", "), 100 * sum(sort(cc, TRUE)[1:5]) / sum(cc), paste(W$location_name[order(-cc)][1:5], collapse = ", "))) }
a <- W %>% arrange(-h_n_2023); cat("Heat deaths 2023, ten largest:", paste(sprintf("%s %s", a$location_name[1:10], fmt(a$h_n_2023[1:10], 0)), collapse = "; "), "\n")
cat(sprintf("India: share of global heat deaths %.1f%% (1990), %.1f%% (2023); South Asia super region %.1f%%, %.1f%%\n", 100 * W$h_n_1990[W$location_name == "India"] / sum(W$h_n_1990), 100 * W$h_n_2023[W$location_name == "India"] / sum(W$h_n_2023),
            100 * A$h_n[A$location_name == "South Asia" & A$year == 1990] / g90$h_n, 100 * A$h_n[A$location_name == "South Asia" & A$year == 2023] / g23$h_n))
a <- W %>% arrange(-h_rate_2023); cat("Highest age standardised heat rates 2023:", paste(sprintf("%s %.1f", a$location_name[1:8], a$h_rate_2023[1:8]), collapse = "; "), "\n")
a <- W %>% arrange(-c_rate_2023); cat("Highest age standardised cold rates 2023:", paste(sprintf("%s %.1f", a$location_name[1:8], a$c_rate_2023[1:8]), collapse = "; "), "\n")
sr <- Wc %>% group_by(region = substr(sr_short[super_region_name], 1, 30)) %>% summarise(n = n(), med_share_1990 = round(median(share_1990), 1), med_share_2023 = round(median(share_2023), 1), heat_dom_1990 = sum(dom90), heat_dom_2023 = sum(dom23), .groups = "drop")
cat("Median country heat share and number of heat dominated countries by super region:\n"); print(as.data.frame(sr), row.names = FALSE)

# ==== 11. Country table ===============================================================================
ct <- W %>% transmute(location_id, location_name, region = sr_short[super_region_name], heat_n_1990 = h_n_1990, heat_n_2023 = h_n_2023, cold_n_1990 = c_n_1990, cold_n_2023 = c_n_2023,
                      heat_rate_1990 = h_rate_1990, heat_rate_2023 = h_rate_2023, cold_rate_1990 = c_rate_1990, cold_rate_2023 = c_rate_2023, heat_frac_1990 = 100 * hp_1990, heat_frac_2023 = 100 * hp_2023, cold_frac_1990 = 100 * cp_1990, cold_frac_2023 = 100 * cp_2023) %>%
  left_join(cmc %>% select(location_id, share90, lo90, hi90, share23, lo23, hi23, p_rise_shared, p_heat_gt_cold23), by = "location_id") %>% left_join(temp %>% filter(year %in% 2016:2020) %>% group_by(location_id) %>% summarise(temp_2016_20 = mean(temp)), by = "location_id") %>%
  left_join(ep %>% select(location_id, y0, y1), by = "location_id")
write_csv(ct, "output/country_table.csv")

# ==== 12. Sensitivity: the same headline results on the series as published ============================
cat("\n==== 12. Sensitivity: published series (cause 729 included, no exclusions), 1990 and 2023 ====\n")
q <- Sp %>% filter(year %in% c(1990, 2023)) %>% mutate(share = 100 * h_rate / (h_rate + c_rate)) %>% select(location_id, year, h_rate, c_rate, share) %>% pivot_wider(names_from = year, values_from = c(h_rate, c_rate, share))
cat(sprintf("Global heat share of counts: %.1f%% to %.1f%% (adjusted %.1f%% to %.1f%%)\n", p90$share_n, p23$share_n, g90$share_n, g23$share_n))
cat(sprintf("Countries where heat exceeds cold: %d to %d (adjusted %d to %d); crossed to heat dominated %d, reverse %d (adjusted %d, %d)\n", sum(q$h_rate_1990 > q$c_rate_1990), sum(q$h_rate_2023 > q$c_rate_2023), sum(Wc$dom90), sum(Wc$dom23),
            sum(q$h_rate_2023 > q$c_rate_2023 & q$h_rate_1990 <= q$c_rate_1990), sum(q$h_rate_2023 <= q$c_rate_2023 & q$h_rate_1990 > q$c_rate_1990), sum(Wc$dom23 & !Wc$dom90), sum(!Wc$dom23 & Wc$dom90)))
cat(sprintf("Median country heat share: %.1f%% to %.1f%% (adjusted %.1f%% to %.1f%%)\n", median(q$share_1990), median(q$share_2023), median(Wc$share_1990), median(Wc$share_2023)))
qb <- Sp %>% mutate(pos = h_rate > 0 & c_rate > 0) %>% group_by(location_id) %>% filter(all(pos)) %>% ungroup() %>% mutate(z = log(h_rate / c_rate)); qs <- qb %>% group_by(year) %>% summarise(sd_z = sd(z))
qw <- qb %>% filter(year %in% c(1990, 2023)) %>% select(location_id, year, z) %>% pivot_wider(names_from = year, values_from = z, names_prefix = "z"); qf <- summary(lm(I(z2023 - z1990) ~ z1990, qw))$coefficients
cat(sprintf("Published series, %d countries with positive rates in all years: SD of log ratio %.3f (1990), %.3f (2023), range over years %.3f to %.3f (adjusted range %.3f to %.3f); beta slope %.4f (SE %.4f); median change in log ratio %+.3f\n", n_distinct(qb$location_id), qs$sd_z[1], qs$sd_z[34], min(qs$sd_z), max(qs$sd_z), min(sig$sd_z), max(sig$sd_z), qf[2, 1], qf[2, 2], median(qw$z2023 - qw$z1990)))
dfc <- Wc %>% select(location_id, location_name, share_1990, share_2023) %>% inner_join(q %>% select(location_id, sp90 = share_1990, sp23 = share_2023), by = "location_id") %>% mutate(d90 = share_1990 - sp90, d23 = share_2023 - sp23)
cat(sprintf("Country heat shares, adjusted minus published: absolute difference above 1 percentage point in %d countries (1990), %d (2023); largest: %s\n", sum(abs(dfc$d90) > 1, na.rm = TRUE), sum(abs(dfc$d23) > 1, na.rm = TRUE),
            paste(sprintf("%s %d (%+.1f)", c(dfc$location_name[order(-abs(dfc$d90))][1:3], dfc$location_name[order(-abs(dfc$d23))][1:4]), rep(c(1990, 2023), c(3, 4)), c(dfc$d90[order(-abs(dfc$d90))][1:3], dfc$d23[order(-abs(dfc$d23))][1:4])), collapse = "; ")))
sink()
cat("analysis complete\n")
