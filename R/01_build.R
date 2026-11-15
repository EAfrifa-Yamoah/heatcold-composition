# 01_build.R
# Assemble the analysis panels from the GBD 2023 extracts in data/raw/ (see R/00_harvest.R).
# Run from the repository root:  Rscript R/01_build.R
#
# Two families of panels are written to data/panel/:
#   (a) "published": heat and cold attributable deaths exactly as GBD 2023 reports them.
#   (b) "adjusted" (suffix _adj): the same quantities after removing the outcome "exposure to forces
#       of nature" (GBD cause 729: earthquakes, storms, floods and other natural disasters), which GBD
#       models as an outcome of high and low temperature. The adjustment subtracts the cause 729
#       component from the heat attributable, cold attributable and all cause series (counts and age
#       standardised rates are additive across causes). Uncertainty bounds are shifted by the same
#       amount, which treats the removed component as known; this is adequate where the component is
#       small and is not relied on in shock years.
# Shock years: country years in which cause 729 accounts for at least 5% of all deaths. In these the
# partition of the remaining causes between heat and cold is distorted (see the manuscript), so they
# are flagged (shock = TRUE) and handled explicitly in R/02_analysis.R.
#
# Outputs:
#   series.csv, series_adj.csv   country x year x sex: heat and cold attributable deaths (count, age
#                                standardised rate, attributable fraction, 95% UI), all cause deaths
#   disaster.csv                 country x year x sex: the cause 729 components
#   shock_years.csv              flagged country years; endpoints.csv: analysis years per country
#   parent.csv, parent_adj.csv   country x {1990, 2023}: non-optimal temperature (parent risk) with UI
#   yll.csv, yll_adj.csv         country x {1990, 2023}: heat and cold attributable YLLs
#   causes.csv, causes_adj.csv   country x {1990, 2023} x 7 level 2 causes (all ages, both sexes)
#   agg_series.csv, agg_series_adj.csv   Global and 7 super regions x year
#   byage.csv, byage_adj.csv     Global and 7 super regions x {1990, 2023} x 25 age bands
#   temperature.csv              country x year population weighted mean temperature (covariate 71)
#   countries.csv                the 204 countries and territories
# Extracts tagged "jsonrounded" hold values rounded to GBD publication precision (three or more
# significant figures); all others are unrounded.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr)})
raw <- "data/raw"; dir.create("data/panel", showWarnings = FALSE)
rd <- function(f) read_csv(file.path(raw, f), show_col_types = FALSE, guess_max = 5000)
rdp <- function(pat) { f <- list.files(raw, pat); stopifnot(length(f) >= 1); map_dfr(f, rd) }
SHOCK_SHARE <- 5   # % of all deaths due to cause 729 that defines a shock year

# long format by sex for one extract: value, lower, upper for rate, count, percent
long_sex <- function(x, id = c("location_id", "location_name")) {
  x %>% select(all_of(id), matches("^(rate|count|percent)_(upper_|lower_)?(male|female|combined)$")) %>%
    pivot_longer(-all_of(id), names_to = "v", values_to = "val") %>%
    mutate(sex = str_extract(v, "(male|female|combined)$"),
           stat = case_when(str_detect(v, "_upper_") ~ "hi", str_detect(v, "_lower_") ~ "lo", TRUE ~ "pt"),
           what = str_extract(v, "^(rate|count|percent)")) %>%
    select(-v) %>% pivot_wider(names_from = c(what, stat), values_from = val)
}

# ---- published country series ---------------------------------------------------------------
years <- 1990:2023
get_series <- function(prefix, tag) map_dfr(years, function(y) {
  x <- rd(sprintf("%s__deaths__as__%d__n204.csv", prefix, y)); stopifnot(nrow(x) == 204)
  long_sex(x) %>% mutate(year = y)
}) %>% rename_with(~ paste0(tag, "_", .x), matches("^(rate|count|percent)_"))
heat <- get_series("top_locations__r337__c294", "h")
cold <- get_series("top_locations__r338__c294", "c")
allc <- get_series("top_locations__c294", "m") %>% select(location_id, year, sex, m_rate = m_rate_pt, d_n = m_count_pt)
series <- heat %>% inner_join(cold %>% select(-location_name), by = c("location_id", "year", "sex")) %>%
  inner_join(allc, by = c("location_id", "year", "sex"))
names(series) <- names(series) %>% str_replace("_rate_pt$", "_rate") %>% str_replace("_count_pt$", "_n") %>% str_replace("_percent_pt$", "_paf") %>%
  str_replace("_count_", "_n_") %>% str_replace("_percent_", "_paf_")
countries <- rd("top_locations__r337__c294__deaths__as__2023__n204.csv") %>% distinct(location_id, location_name, super_region_name) %>% arrange(super_region_name, location_name)
stopifnot(nrow(countries) == 204, nrow(series) == 204 * 34 * 3, !anyNA(series$h_rate), !anyNA(series$c_rate), !anyNA(series$m_rate))
series <- series %>% left_join(countries %>% select(location_id, super_region_name), by = "location_id") %>%
  select(location_id, location_name, super_region_name, year, sex, everything())
write_csv(series, "data/panel/series.csv"); write_csv(countries, "data/panel/countries.csv")

# ---- cause 729 components (yearly, all locations) ---------------------------------------------
fon <- function(pat, tag) rdp(pat) %>% distinct(location_id, year, .keep_all = TRUE) %>%
  select(location_id, year, matches("^(count|rate_age_std|rate|percent)_(male|female|combined)$")) %>%
  pivot_longer(-c(location_id, year), names_to = c("what", "sex"), names_pattern = "^(count|rate_age_std|rate|percent)_(male|female|combined)$") %>%
  pivot_wider(names_from = what, values_from = value) %>%
  transmute(location_id, year, sex, n = count, asr = rate_age_std, crude = rate, paf = percent) %>% rename_with(~ paste0(tag, "_", .x), c(n, asr, crude, paf))
dis_all <- fon("^multi_location_trend__r337__c729__", "h729") %>% full_join(fon("^multi_location_trend__r338__c729__", "c729"), by = c("location_id", "year", "sex")) %>%
  full_join(fon("^multi_location_trend__c729__", "d729") %>% select(-d729_paf), by = c("location_id", "year", "sex"))
agg_ids <- rd("multi_location_trend__c294__deaths__all__1990to2023__n8.csv") %>% distinct(location_id) %>% pull()
dis <- dis_all %>% filter(location_id %in% countries$location_id); agg_dis <- dis_all %>% filter(location_id %in% agg_ids, sex == "combined")
stopifnot(nrow(dis) == 204 * 34 * 3, nrow(agg_dis) == 8 * 34, !anyNA(dis$h729_n), !anyNA(dis$c729_n), !anyNA(dis$d729_n), !anyNA(dis$h729_asr), !anyNA(dis$c729_asr), !anyNA(dis$d729_asr))
write_csv(dis, "data/panel/disaster.csv")

# ---- adjusted country series and shock years ----------------------------------------------------
series_adj <- series %>% inner_join(dis, by = c("location_id", "year", "sex")) %>%
  mutate(across(c(h_rate, h_rate_hi, h_rate_lo), ~ .x - h729_asr), across(c(h_n, h_n_hi, h_n_lo), ~ .x - h729_n),
         across(c(c_rate, c_rate_hi, c_rate_lo), ~ .x - c729_asr), across(c(c_n, c_n_hi, c_n_lo), ~ .x - c729_n),
         d729_share = 100 * d729_n / d_n, m_rate = m_rate - d729_asr, d_n = d_n - d729_n,
         h_paf = 100 * h_rate / m_rate, c_paf = 100 * c_rate / m_rate) %>%
  select(-matches("_paf_(hi|lo)$"), -matches("729_(n|asr|crude|paf)$"))
shock <- series_adj %>% filter(sex == "combined") %>% transmute(location_id, year, shock = d729_share >= SHOCK_SHARE)
series_adj <- series_adj %>% left_join(shock, by = c("location_id", "year"))
stopifnot(nrow(series_adj) == nrow(series), all(series_adj$d_n > 0), all(series_adj$m_rate > 0))
write_csv(series_adj, "data/panel/series_adj.csv")
sy <- series %>% filter(sex == "combined") %>% inner_join(dis %>% filter(sex == "combined"), by = c("location_id", "year")) %>% inner_join(shock, by = c("location_id", "year")) %>% filter(shock) %>%
  transmute(location_id, location_name, year, d729_n, d729_share = 100 * d729_n / d_n, heat_published = h_n, heat_729 = h729_n, heat_729_paf = h729_paf, cold_published = c_n, cold_729 = c729_n, cold_729_paf = c729_paf) %>% arrange(desc(d729_share))
write_csv(sy, "data/panel/shock_years.csv")
# analysis years at the two ends of the series: 1990 and 2023, or the nearest year that is not a shock year
ok <- shock %>% filter(!shock)
endpoints <- countries %>% select(location_id, location_name) %>% rowwise() %>%
  mutate(y0 = min(ok$year[ok$location_id == location_id]), y1 = max(ok$year[ok$location_id == location_id])) %>% ungroup()
write_csv(endpoints, "data/panel/endpoints.csv")

# ---- parent risk and YLLs (1990 and 2023) ---------------------------------------------------------
parent <- map_dfr(c(1990, 2023), function(y) long_sex(rd(sprintf("top_locations__r331__c294__deaths__as__%d__n204.csv", y))) %>% mutate(year = y))
names(parent) <- names(parent) %>% str_replace("^(rate|count|percent)_", "p_\\1_")
write_csv(parent, "data/panel/parent.csv")
parent_adj <- parent %>% inner_join(dis, by = c("location_id", "year", "sex")) %>% mutate(across(c(p_rate_pt, p_rate_hi, p_rate_lo), ~ .x - h729_asr - c729_asr)) %>%
  select(location_id, location_name, sex, year, p_rate_pt, p_rate_hi, p_rate_lo)
write_csv(parent_adj, "data/panel/parent_adj.csv")
yll <- map_dfr(c(1990, 2023), function(y) {
  a <- long_sex(rd(sprintf("top_locations__r337__c294__ylls__as__%d__n204.csv", y))) %>% transmute(location_id, sex, h_yll = count_pt, h_yll_rate = rate_pt)
  b <- long_sex(rd(sprintf("top_locations__r338__c294__ylls__as__%d__n204.csv", y))) %>% transmute(location_id, sex, c_yll = count_pt, c_yll_rate = rate_pt)
  inner_join(a, b, by = c("location_id", "sex")) %>% mutate(year = y)
})
write_csv(yll, "data/panel/yll.csv")
yll729 <- map_dfr(c(1990, 2023), function(y) {
  a <- rd(sprintf("top_locations__r337__c729__ylls__as__%d__n204__jsonrounded.csv", y)) %>% transmute(location_id, h729_yll = count_combined)
  b <- rd(sprintf("top_locations__r338__c729__ylls__as__%d__n204__jsonrounded.csv", y)) %>% transmute(location_id, c729_yll = count_combined)
  full_join(a, b, by = "location_id") %>% mutate(year = y)
})
yll_adj <- yll %>% filter(sex == "combined") %>% left_join(yll729, by = c("location_id", "year")) %>%
  mutate(h_yll = h_yll - coalesce(h729_yll, 0), c_yll = c_yll - coalesce(c729_yll, 0)) %>% select(location_id, year, h_yll, c_yll)
stopifnot(nrow(yll_adj) == 408)
write_csv(yll_adj, "data/panel/yll_adj.csv")

# ---- cause pairs (level 2 causes, all ages, both sexes) ---------------------------------------
cause_lab <- c(`491` = "Cardiovascular diseases", `508` = "Chronic respiratory diseases", `974` = "Diabetes and kidney diseases",
               `956` = "Respiratory infections and tuberculosis", `688` = "Transport injuries", `696` = "Unintentional injuries",
               `717` = "Self-harm and interpersonal violence")
causes <- map_dfr(names(cause_lab), function(cid) map_dfr(c(1990, 2023), function(y) {
  h <- rd(sprintf("top_locations__r337__c%s__deaths__all__%d__n204.csv", cid, y)) %>% transmute(location_id, h_n = count_combined, h_paf = percent_combined)
  k <- rd(sprintf("top_locations__r338__c%s__deaths__all__%d__n204.csv", cid, y)) %>% transmute(location_id, c_n = count_combined, c_paf = percent_combined)
  d <- rd(sprintf("top_locations__c%s__deaths__all__%d__n204.csv", cid, y)) %>% transmute(location_id, d_n = count_combined)
  stopifnot(nrow(h) == 204, nrow(k) == 204, nrow(d) == 204)
  h %>% inner_join(k, by = "location_id") %>% inner_join(d, by = "location_id") %>% mutate(year = y, cause_id = as.integer(cid), cause = cause_lab[[cid]])
}))
stopifnot(nrow(causes) == 204 * 2 * 7)
write_csv(causes, "data/panel/causes.csv")
causes_adj <- causes %>% left_join(dis %>% filter(sex == "combined") %>% select(location_id, year, h729_n, c729_n, d729_n), by = c("location_id", "year")) %>%
  mutate(u = cause_id == 696, h_n = h_n - ifelse(u, h729_n, 0), c_n = c_n - ifelse(u, c729_n, 0), d_n = d_n - ifelse(u, d729_n, 0),
         h_paf = 100 * h_n / d_n, c_paf = 100 * c_n / d_n, cause = ifelse(u, "Unintentional injuries excluding forces of nature", cause)) %>% select(-u, -h729_n, -c729_n, -d729_n)
write_csv(causes_adj, "data/panel/causes_adj.csv")

# ---- aggregates: yearly series and age patterns ------------------------------------------------
ag1 <- function(f, tag) rd(f) %>% transmute(location_id, location_name, year, !!paste0(tag, "_n") := count_combined, !!paste0(tag, "_n_lo") := count_lower_combined,
                                           !!paste0(tag, "_n_hi") := count_upper_combined, !!paste0(tag, "_crude") := rate_combined, !!paste0(tag, "_asr") := rate_age_std_combined,
                                           !!paste0(tag, "_paf") := percent_combined, !!paste0(tag, "_paf_lo") := percent_lower_combined, !!paste0(tag, "_paf_hi") := percent_upper_combined)
agg <- ag1("multi_location_trend__r337__c294__deaths__all__1990to2023__n8.csv", "h") %>%
  inner_join(ag1("multi_location_trend__r338__c294__deaths__all__1990to2023__n8.csv", "c") %>% select(-location_name), by = c("location_id", "year")) %>%
  inner_join(ag1("multi_location_trend__r331__c294__deaths__all__1990to2023__n8.csv", "p") %>% select(-location_name), by = c("location_id", "year")) %>%
  inner_join(rd("multi_location_trend__c294__deaths__all__1990to2023__n8.csv") %>% transmute(location_id, year, d_n = count_combined, m_crude = rate_combined, m_asr = rate_age_std_combined), by = c("location_id", "year"))
stopifnot(nrow(agg) == 8 * 34)
write_csv(agg, "data/panel/agg_series.csv")
agg_adj <- agg %>% inner_join(agg_dis %>% select(-sex), by = c("location_id", "year")) %>%
  mutate(across(c(h_n, h_n_lo, h_n_hi), ~ .x - h729_n), across(c(c_n, c_n_lo, c_n_hi), ~ .x - c729_n), across(c(p_n, p_n_lo, p_n_hi), ~ .x - h729_n - c729_n),
         h_crude = h_crude - h729_crude, c_crude = c_crude - c729_crude, p_crude = p_crude - h729_crude - c729_crude,
         h_asr = h_asr - h729_asr, c_asr = c_asr - c729_asr, p_asr = p_asr - h729_asr - c729_asr,
         d_n = d_n - d729_n, m_crude = m_crude - d729_crude, m_asr = m_asr - d729_asr,
         across(c(h_n, h_n_lo, h_n_hi, c_n, c_n_lo, c_n_hi, p_n, p_n_lo, p_n_hi), ~ 100 * .x / d_n, .names = "{.col}_f")) %>%
  mutate(h_paf = h_n_f, h_paf_lo = h_n_lo_f, h_paf_hi = h_n_hi_f, c_paf = c_n_f, c_paf_lo = c_n_lo_f, c_paf_hi = c_n_hi_f, p_paf = p_n_f, p_paf_lo = p_n_lo_f, p_paf_hi = p_n_hi_f) %>%
  select(-ends_with("_f"), -matches("729_(asr|crude|paf)$")) %>% rename(h729_n_removed = h729_n, c729_n_removed = c729_n, d729_n_removed = d729_n)
write_csv(agg_adj, "data/panel/agg_series_adj.csv")
byage <- map_dfr(c(1990, 2023), function(y) {
  h <- rd(sprintf("risk_attributable_burden__high_temperature__all_causes__deaths__byage__%d__n8.csv", y)) %>% transmute(location_id, location_name, age_group_id, age_group, h_n = count_combined, h_paf = percent_combined)
  k <- rd(sprintf("risk_attributable_burden__low_temperature__all_causes__deaths__byage__%d__n8.csv", y)) %>% transmute(location_id, age_group_id, c_n = count_combined, c_paf = percent_combined)
  d <- rd(sprintf("cause_burden__call_causes__deaths__byage__%d__n8.csv", y)) %>% transmute(location_id, age_group_id, d_n = count_combined, d_rate = rate_combined)
  h %>% inner_join(k, by = c("location_id", "age_group_id")) %>% inner_join(d, by = c("location_id", "age_group_id")) %>% mutate(year = y, pop = d_n / d_rate * 1e5)
})
stopifnot(nrow(byage) == 8 * 25 * 2)
write_csv(byage, "data/panel/byage.csv")
byage729 <- map_dfr(c(1990, 2023), function(y) {
  h <- rd(sprintf("risk_attributable_burden__high_temperature__exposure_to_forces_of_nature__deaths__byage__%d__n8__jsonrounded.csv", y)) %>% transmute(location_id, age_group_id, h729 = count_combined)
  k <- rd(sprintf("risk_attributable_burden__low_temperature__exposure_to_forces_of_nature__deaths__byage__%d__n8__jsonrounded.csv", y)) %>% transmute(location_id, age_group_id, c729 = count_combined)
  d <- rd(sprintf("cause_burden__cexposure_to_forces_of_nature__deaths__byage__%d__n8__jsonrounded.csv", y)) %>% transmute(location_id, age_group_id, d729 = count_combined)
  full_join(h, k, by = c("location_id", "age_group_id")) %>% full_join(d, by = c("location_id", "age_group_id")) %>% mutate(year = y)
})
byage_adj <- byage %>% left_join(byage729, by = c("location_id", "age_group_id", "year")) %>%
  mutate(h_n = h_n - coalesce(h729, 0), c_n = c_n - coalesce(c729, 0), d_n = d_n - coalesce(d729, 0), h_paf = 100 * h_n / d_n, c_paf = 100 * c_n / d_n) %>% select(-h729, -c729, -d729, -d_rate)
write_csv(byage_adj, "data/panel/byage_adj.csv")

# ---- temperature covariate ---------------------------------------------------------------------
temp <- map_dfr(list.files(raw, "^cov__population_weighted_mean_temperature", full.names = TRUE), ~ read_csv(.x, show_col_types = FALSE)) %>%
  transmute(location_id, year, temp = value_combined) %>% distinct()
stopifnot(nrow(temp) == 204 * 34, setequal(temp$location_id, countries$location_id))
# the covariate is carried forward unchanged after 2020; flag those years
frozen <- temp %>% arrange(location_id, year) %>% group_by(location_id) %>% mutate(same = temp == lag(temp)) %>% group_by(year) %>% summarise(share_same = mean(same, na.rm = TRUE))
temp <- temp %>% mutate(carried_forward = year >= 2021)
write_csv(temp, "data/panel/temperature.csv")

# ---- build report ----------------------------------------------------------------------------------
sink("data/panel/_build_report.txt")
cat("series rows:", nrow(series), " causes rows:", nrow(causes), " agg rows:", nrow(agg), " byage rows:", nrow(byage), " temperature rows:", nrow(temp), "\n")
cat("share of countries whose temperature equals the previous year, by year (last 6):\n"); print(as.data.frame(tail(frozen, 6)))
chk <- series %>% filter(sex == "combined", year %in% c(1990, 2023)) %>% group_by(year) %>% summarise(heat = sum(h_n), cold = sum(c_n), all = sum(d_n))
g <- agg %>% filter(location_id == 1, year %in% c(1990, 2023)) %>% select(year, h_n, c_n, d_n)
cat("\nsum over 204 countries vs published Global (deaths), relative difference:\n"); print(as.data.frame(left_join(chk, g, by = "year") %>% transmute(year, heat_rel = heat / h_n - 1, cold_rel = cold / c_n - 1, all_rel = all / d_n - 1)))
c2 <- dis %>% filter(sex == "combined") %>% group_by(year) %>% summarise(h = sum(h729_n), c = sum(c729_n), d = sum(d729_n)) %>% inner_join(agg_dis %>% filter(location_id == 1) %>% select(year, h729_n, c729_n, d729_n), by = "year")
cat("\ncause 729: sum over 204 countries vs Global, range of relative differences over 34 years: heat", paste(signif(range(c2$h / c2$h729_n - 1), 2), collapse = " to "), "; cold", paste(signif(range(c2$c / c2$c729_n - 1), 2), collapse = " to "),
    "; deaths", paste(signif(range(c2$d / c2$d729_n - 1), 2), collapse = " to "), "\n")
c3 <- byage729 %>% group_by(location_id, year) %>% summarise(h = sum(h729, na.rm = TRUE), c = sum(c729, na.rm = TRUE), d = sum(d729, na.rm = TRUE), .groups = "drop") %>% inner_join(agg_dis %>% select(location_id, year, h729_n, c729_n, d729_n), by = c("location_id", "year"))
cat("cause 729 by age: sum over ages vs all ages, max absolute relative difference: heat", signif(max(abs(c3$h / c3$h729_n - 1)), 2), "; cold", signif(max(abs(c3$c / c3$c729_n - 1)), 2), "; deaths", signif(max(abs(c3$d / c3$d729_n - 1)), 2), "\n")
cat("\nshock years (cause 729 at least", SHOCK_SHARE, "% of all deaths):", nrow(sy), "country years in", n_distinct(sy$location_id), "countries\n"); print(as.data.frame(sy %>% mutate(across(where(is.numeric) & !c(location_id, year), ~ round(.x, 1)))), row.names = FALSE)
cat("\nendpoints other than 1990 and 2023:\n"); print(as.data.frame(endpoints %>% filter(y0 != 1990 | y1 != 2023)), row.names = FALSE)
cat("\nadjusted series: country years with non-positive heat rate", sum(series_adj$h_rate[series_adj$sex == "combined"] <= 0), ", cold rate", sum(series_adj$c_rate[series_adj$sex == "combined"] <= 0),
    "; published: heat", sum(series$h_rate[series$sex == "combined"] <= 0), ", cold", sum(series$c_rate[series$sex == "combined"] <= 0), "\n")
sink()
cat(readLines("data/panel/_build_report.txt"), sep = "\n")
