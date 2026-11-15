# 99_independent_check.R
# Recomputes headline results directly from the extracts in data/raw/ with base R only, by a code
# path that shares nothing with R/01_build.R and R/02_analysis.R (no panels, no dplyr, Shapley values
# by enumerating orderings rather than subsets), and compares them with the values written by the
# main pipeline. Stops with an error if any comparison fails. Run after R/02_analysis.R.
raw <- "data/raw"; rd <- function(f) read.csv(file.path(raw, f), check.names = FALSE, encoding = "UTF-8")
chk <- function(label, mine, theirs, tol = 1e-6) { ok <- all(abs(mine - theirs) <= tol * pmax(1, abs(theirs))); cat(sprintf("%-78s %s  (check %s; pipeline %s)\n", label, ifelse(ok, "OK  ", "FAIL"), paste(signif(mine, 7), collapse = ", "), paste(signif(theirs, 7), collapse = ", "))); if (!ok) stop("mismatch: ", label) }
bindf <- function(pat) { fs <- list.files(raw, pat, full.names = TRUE); x <- do.call(rbind, lapply(fs, function(f) read.csv(f, check.names = FALSE)[, c("location_id", "year", "count_combined", "rate_age_std_combined", "percent_combined")])); x[!duplicated(x[, c("location_id", "year")]), ] }

# 1. Global totals, adjusted for cause 729, from the aggregate extracts ------------------------------------
gl <- function(f) { x <- rd(f); x[x$location_id == 1, ] }
H <- gl("multi_location_trend__r337__c294__deaths__all__1990to2023__n8.csv"); K <- gl("multi_location_trend__r338__c294__deaths__all__1990to2023__n8.csv")
h7 <- bindf("^multi_location_trend__r337__c729__"); c7 <- bindf("^multi_location_trend__r338__c729__"); d7 <- bindf("^multi_location_trend__c729__")
g <- function(x, y) x$count_combined[x$location_id == 1 & x$year == y]
adj <- sapply(c(1990, 2023), function(y) c(heat = g(H, y) - g(h7, y), cold = g(K, y) - g(c7, y)))
pipe <- read.csv("output/table1_regions.csv"); pg <- pipe[pipe$location_id == 1, ]
chk("Global heat deaths excluding forces of nature, 1990 and 2023", adj["heat", ], pg$h_n[order(pg$year)])
chk("Global cold deaths excluding forces of nature, 1990 and 2023", adj["cold", ], pg$c_n[order(pg$year)])
chk("Global heat share (%), 1990 and 2023", 100 * adj["heat", ] / colSums(adj), pg$share_n[order(pg$year)])

# 2. Country series: adjusted rates, shock years, classification -------------------------------------------
yr <- 1990:2023
ser <- do.call(rbind, lapply(yr, function(y) { a <- rd(sprintf("top_locations__r337__c294__deaths__as__%d__n204.csv", y)); b <- rd(sprintf("top_locations__r338__c294__deaths__as__%d__n204.csv", y)); m <- rd(sprintf("top_locations__c294__deaths__as__%d__n204.csv", y))
  b <- b[match(a$location_id, b$location_id), ]; m <- m[match(a$location_id, m$location_id), ]
  data.frame(id = a$location_id, year = y, h = a$rate_combined, cc = b$rate_combined, m = m$rate_combined, d = m$count_combined) }))
key <- function(x) paste(x$location_id, x$year); k <- paste(ser$id, ser$year)
ser$h <- ser$h - h7$rate_age_std_combined[match(k, key(h7))]; ser$cc <- ser$cc - c7$rate_age_std_combined[match(k, key(c7))]
ser$share729 <- 100 * d7$count_combined[match(k, key(d7))] / ser$d; ser$shock <- ser$share729 >= 5
ser$m <- ser$m - d7$rate_age_std_combined[match(k, key(d7))]   # all cause rate excluding cause 729
stopifnot(!anyNA(ser$h), !anyNA(ser$cc), nrow(ser) == 204 * 34)
chk("Number of shock country years and countries", c(sum(ser$shock), length(unique(ser$id[ser$shock]))), c(nrow(read.csv("data/panel/shock_years.csv")), length(unique(read.csv("data/panel/shock_years.csv")$location_id))))
endv <- function(first) do.call(rbind, lapply(split(ser[!ser$shock, ], ser$id[!ser$shock]), function(x) x[if (first) which.min(x$year) else which.max(x$year), ]))
e0 <- endv(TRUE); e1 <- endv(FALSE); stopifnot(identical(e0$id, e1$id), nrow(e0) == 204)
ct <- read.csv("output/country_table.csv"); ct <- ct[match(e0$id, ct$location_id), ]
chk("Countries where heat exceeds cold, 1990 and 2023", c(sum(e0$h > e0$cc), sum(e1$h > e1$cc)), c(sum(ct$heat_rate_1990 > ct$cold_rate_1990), sum(ct$heat_rate_2023 > ct$cold_rate_2023)))
s0 <- 100 * e0$h / (e0$h + e0$cc); s1 <- 100 * e1$h / (e1$h + e1$cc)
chk("Median country heat share (%), 1990 and 2023", c(median(s0), median(s1)), c(median(ct$share90), median(ct$share23)))
chk("Country heat rates in 2023 (max absolute difference from the pipeline)", max(abs(e1$h - ct$heat_rate_2023)), 0, tol = 1e-8)
pos <- e0$h > 0 & e0$cc > 0 & e1$h > 0 & e1$cc > 0; ci <- read.csv("output/country_identity.csv")
chk("Countries with positive rates at both ends; median change in log heat and cold fractions", c(sum(pos), median(log((e1$h / e1$m) / (e0$h / e0$m))[pos]), median(log((e1$cc / e1$m) / (e0$cc / e0$m))[pos])), c(nrow(ci), median(ci$dHp), median(ci$dCp)))
# dispersion of the log ratio in the balanced set (shock years interpolated)
okc <- sapply(split(ser, ser$id), function(x) all((x$h > 0 & x$cc > 0)[!x$shock])); bal <- as.integer(names(okc)[okc])
zs <- sapply(bal, function(i) { x <- ser[ser$id == i, ]; x <- x[order(x$year), ]; z <- suppressWarnings(log(x$h / x$cc)); approx(x$year[!x$shock], z[!x$shock], xout = yr, rule = 2)$y })
sg <- read.csv("output/sigma_series.csv")
chk("Balanced set size; SD of log ratio in 1990 and 2023", c(length(bal), sd(zs[1, ]), sd(zs[34, ])), c(193, sg$sd_z[sg$year == 1990], sg$sd_z[sg$year == 2023]))
cat(sprintf("%-78s %s\n", "Median change in log ratio 1990 to 2023 (balanced set; compare stats.txt section 4)", signif(median(zs[34, ] - zs[1, ]), 4)))

# 3. Shapley decomposition of the global heat count by enumerating the 24 orderings ------------------------------
age <- function(f, col = "count_combined") { x <- rd(f); x <- x[x$location_id == 1, ]; x[order(x$age_group_id), c("age_group_id", col)] }
mk <- function(y) { h <- age(sprintf("risk_attributable_burden__high_temperature__all_causes__deaths__byage__%d__n8.csv", y)); d <- age(sprintf("cause_burden__call_causes__deaths__byage__%d__n8.csv", y)); r <- age(sprintf("cause_burden__call_causes__deaths__byage__%d__n8.csv", y), "rate_combined")
  h729 <- age(sprintf("risk_attributable_burden__high_temperature__exposure_to_forces_of_nature__deaths__byage__%d__n8__jsonrounded.csv", y)); d729 <- age(sprintf("cause_burden__cexposure_to_forces_of_nature__deaths__byage__%d__n8__jsonrounded.csv", y))
  stopifnot(identical(h$age_group_id, d$age_group_id)); pop <- d$count_combined / r$rate_combined * 1e5
  hh <- h$count_combined - ifelse(is.na(m <- match(h$age_group_id, h729$age_group_id)), 0, h729$count_combined[m]); dd <- d$count_combined - ifelse(is.na(m <- match(d$age_group_id, d729$age_group_id)), 0, d729$count_combined[m])
  list(N = sum(pop), s = pop / sum(pop), m = dd / pop, p = hh / dd) }
x0 <- mk(1990); x1 <- mk(2023); f <- function(x) x$N * sum(x$s * x$m * x$p)
perms <- function(v) if (length(v) == 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
phi <- setNames(numeric(4), names(x0)); for (p in perms(names(x0))) { x <- x0; for (nm in p) { before <- f(x); x[[nm]] <- x1[[nm]]; phi[nm] <- phi[nm] + (f(x) - before) / 24 } }
dc <- read.csv("output/decomposition_counts.csv"); dg <- dc[dc$location_id == 1 & dc$risk == "h", ]
chk("Shapley decomposition of global heat deaths (population, age, mortality, fraction)", unname(phi), c(dg$pop, dg$age, dg$mort, dg$paf), tol = 1e-6)

# 4. Cause decomposition of the global heat fraction ---------------------------------------------------------------
cid <- c(491, 508, 974, 956, 688, 696, 717)
cs <- function(y) t(sapply(cid, function(cc) { h <- rd(sprintf("top_locations__r337__c%d__deaths__all__%d__n204.csv", cc, y)); d <- rd(sprintf("top_locations__c%d__deaths__all__%d__n204.csv", cc, y)); c(h = sum(h$count_combined), d = sum(d$count_combined)) }))
tot <- function(y) sum(rd(sprintf("top_locations__c294__deaths__as__%d__n204.csv", y))$count_combined) - sum(d7$count_combined[d7$year == y & d7$location_id %in% ser$id])
adj7 <- function(m, y) { m[cid == 696, "h"] <- m[cid == 696, "h"] - sum(h7$count_combined[h7$year == y & h7$location_id %in% ser$id]); m[cid == 696, "d"] <- m[cid == 696, "d"] - sum(d7$count_combined[d7$year == y & d7$location_id %in% ser$id]); m }
a0 <- adj7(cs(1990), 1990); a1 <- adj7(cs(2023), 2023); w0 <- a0[, "d"] / tot(1990); w1 <- a1[, "d"] / tot(2023); p0 <- a0[, "h"] / a0[, "d"]; p1 <- a1[, "h"] / a1[, "d"]
gw <- read.csv("output/global_causes.csv"); gw <- gw[match(c("Cardiovascular diseases", "Chronic respiratory diseases", "Diabetes and kidney diseases", "Respiratory infections and tuberculosis", "Transport injuries", "Unintentional injuries excluding forces of nature", "Self-harm and interpersonal violence"), gw$cause), ]
chk("Global heat fraction: within cause and cause mix components (percentage points)", 100 * c(sum((p1 - p0) * (w0 + w1) / 2), sum((w1 - w0) * (p0 + p1) / 2)), 100 * c(sum((gw$ph_2023 - gw$ph_1990) * (gw$w_1990 + gw$w_2023) / 2), sum((gw$w_2023 - gw$w_1990) * (gw$ph_1990 + gw$ph_2023) / 2)))
cat("\nAll independent checks passed.\n")
