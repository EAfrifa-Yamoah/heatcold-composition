# 03_figures.R
# Figures for the manuscript and supplement, drawn from output/*.csv (run R/02_analysis.R first).
# Conventions: no titles inside figures (captions are in the text); one colour, one meaning across
# the whole document (Okabe-Ito palette); a single shared legend per figure.
#   heat = vermilion, cold = blue; decomposition components: population growth = grey, age structure
#   = orange, all cause mortality = bluish green, attributable fraction = reddish purple.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(sf); library(scales); library(ggrepel)})
dir.create("figures", showWarnings = FALSE); invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
ro <- function(f) read_csv(file.path("output", f), show_col_types = FALSE)
HEAT <- "#D55E00"; COLD <- "#0072B2"; GREY <- "#999999"; ORANGE <- "#E69F00"; GREEN <- "#009E73"; PURPLE <- "#CC79A7"; INK <- "#222222"
risk_cols <- c(Heat = HEAT, Cold = COLD)
theme_set(theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"), legend.position = "bottom", legend.key.height = unit(3.5, "mm"),
                                               plot.tag = element_text(face = "bold", size = 10), axis.title = element_text(size = 8.5), strip.text = element_text(face = "bold", size = 8.5, hjust = 0), plot.margin = margin(4, 6, 4, 4)))
save_fig <- function(p, name, w, h) { ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf); ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white") }
reg_order <- c("Global", "High income", "Central and eastern Europe, central Asia", "Latin America and Caribbean", "North Africa and Middle East", "South Asia", "Southeast and east Asia, Oceania", "Sub-Saharan Africa")

# ---- Figure 1: global series ---------------------------------------------------------------------
A <- ro("agg_shares.csv"); g <- A %>% filter(location_id == 1)
d1 <- bind_rows(g %>% transmute(year, risk = "Heat", series = "Excluding forces of nature", n = h_n), g %>% transmute(year, risk = "Cold", series = "Excluding forces of nature", n = c_n),
                g %>% transmute(year, risk = "Heat", series = "As published", n = h_pub), g %>% transmute(year, risk = "Cold", series = "As published", n = c_pub)) %>% mutate(series = factor(series, c("Excluding forces of nature", "As published")))
p1a <- ggplot(d1, aes(year, n / 1e3, colour = risk, linetype = series)) + geom_line(linewidth = 0.55) + scale_colour_manual(values = risk_cols, name = NULL) + scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_y_continuous(limits = c(0, NA), labels = comma) + labs(x = NULL, y = "Attributable deaths (thousands)", tag = "A")
d1b <- bind_rows(g %>% transmute(year, series = "Excluding forces of nature", share = share_n), g %>% transmute(year, series = "As published", share = share_pub)) %>% mutate(series = factor(series, c("Excluding forces of nature", "As published")))
p1b <- ggplot(d1b, aes(year, share, linetype = series)) + geom_line(linewidth = 0.55, colour = INK) + scale_linetype_manual(values = c("solid", "22"), name = NULL, guide = "none") + labs(x = NULL, y = "Heat share of temperature\nattributable deaths (%)", tag = "B")
d1c <- g %>% transmute(year, Heat = 100 * (h_asr / m_asr) / (h_asr[year == 1990] / m_asr[year == 1990]), Cold = 100 * (c_asr / m_asr) / (c_asr[year == 1990] / m_asr[year == 1990])) %>% pivot_longer(-year, names_to = "risk", values_to = "idx")
p1c <- ggplot(d1c, aes(year, idx, colour = risk)) + geom_hline(yintercept = 100, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.55) + scale_colour_manual(values = risk_cols, guide = "none") +
  labs(x = NULL, y = "Age standardised attributable\nfraction (1990 = 100)", tag = "C")
d1d <- g %>% transmute(year, Heat = 100 * h_asr / h_asr[year == 1990], Cold = 100 * c_asr / c_asr[year == 1990]) %>% pivot_longer(-year, names_to = "risk", values_to = "idx")
p1d <- ggplot(d1d, aes(year, idx, colour = risk)) + geom_hline(yintercept = 100, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.55) + scale_colour_manual(values = risk_cols, guide = "none") +
  labs(x = NULL, y = "Age standardised attributable\ndeath rate (1990 = 100)", tag = "D")
save_fig((p1a + p1b) / (p1c + p1d) + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "fig1_global_series", 7.2, 5.6)

# ---- Figure 2: maps --------------------------------------------------------------------------------
ct <- ro("country_table.csv")
world <- st_as_sf(rnaturalearthdata::countries110) %>% filter(admin != "Antarctica") %>% select(admin, geometry)
recode_ne <- c("Russia" = "Russian Federation", "The Bahamas" = "Bahamas", "East Timor" = "Timor-Leste", "Bolivia" = "Bolivia (Plurinational State of)", "Venezuela" = "Venezuela (Bolivarian Republic of)", "Ivory Coast" = "Côte d'Ivoire",
               "Republic of the Congo" = "Congo", "eSwatini" = "Eswatini", "Laos" = "Lao People's Democratic Republic", "Vietnam" = "Viet Nam", "North Korea" = "Democratic People's Republic of Korea", "South Korea" = "Republic of Korea",
               "Iran" = "Iran (Islamic Republic of)", "Syria" = "Syrian Arab Republic", "Moldova" = "Republic of Moldova", "Turkey" = "Türkiye", "Brunei" = "Brunei Darussalam", "Republic of Serbia" = "Serbia", "Somaliland" = "Somalia", "Northern Cyprus" = "Cyprus", "Kosovo" = "Serbia")
world <- world %>% mutate(location_name = enc2utf8(unname(ifelse(admin %in% names(recode_ne), recode_ne[admin], admin)))) %>% left_join(ct %>% select(location_name, share90, share23), by = "location_name") %>%
  mutate(share23c = pmin(pmax(share23, 0), 100), dshare = pmin(pmax(pmin(pmax(share23, 0), 100) - pmin(pmax(share90, 0), 100), -40), 40)) %>% st_transform("+proj=robin")
cat("map polygons without a GBD match:", paste(world$admin[is.na(world$share23)], collapse = "; "), "\n"); cat("GBD countries not drawn at this map resolution:", length(setdiff(ct$location_name, world$location_name)), "\n")
map_theme <- theme(axis.text = element_blank(), panel.grid.major = element_line(linewidth = 0.15, colour = "grey90"), legend.position = "right", legend.key.height = unit(9, "mm"), legend.key.width = unit(3, "mm"), legend.title = element_text(size = 7.5))
p2a <- ggplot(world) + geom_sf(aes(fill = share23c), colour = "white", linewidth = 0.1) + scale_fill_gradient2(low = COLD, mid = "grey96", high = HEAT, midpoint = 50, limits = c(0, 100), na.value = "grey80", name = "Heat share\n2023 (%)") + labs(tag = "A") + map_theme
p2b <- ggplot(world) + geom_sf(aes(fill = dshare), colour = "white", linewidth = 0.1) + scale_fill_gradient2(low = COLD, mid = "grey96", high = HEAT, midpoint = 0, limits = c(-40, 40), na.value = "grey80", name = "Change since\n1990 (points)", breaks = c(-40, -20, 0, 20, 40), labels = c("-40", "-20", "0", "+20", "+40 or more")) + labs(tag = "B") + map_theme
save_fig(p2a / p2b & theme(plot.margin = margin(0, 2, 0, 2)), "fig2_maps", 7.2, 6.3)

# ---- Figure 3: drivers of attributable deaths and rates ------------------------------------------------
dec <- ro("decomposition_counts.csv") %>% mutate(region = factor(region, rev(reg_order)), risk = factor(ifelse(risk == "h", "Heat", "Cold"), c("Heat", "Cold")))
dl <- dec %>% select(region, risk, pop_pct, age_pct, mort_pct, paf_pct) %>% pivot_longer(-c(region, risk), names_to = "comp", values_to = "v") %>%
  mutate(comp = factor(comp, c("pop_pct", "age_pct", "mort_pct", "paf_pct"), c("Population growth", "Age structure", "Age specific all cause mortality", "Age specific attributable fraction")))
p3a <- ggplot(dl, aes(v, region, fill = comp)) + geom_col(width = 0.68) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(data = dec, aes(net_pct, region), inherit.aes = FALSE, shape = 23, fill = "white", colour = INK, size = 1.9, stroke = 0.6) +
  facet_wrap(~ risk, scales = "free_x") + scale_fill_manual(values = c(GREY, ORANGE, GREEN, PURPLE), name = NULL) + labs(x = "Contribution to the change in attributable deaths, 1990 to 2023 (% of the 1990 count)", y = NULL, tag = "A") +
  guides(fill = guide_legend(nrow = 2))
P <- ro("country_identity.csv")
pl <- bind_rows(P %>% transmute(risk = "Heat", dM, dF = dHp), P %>% transmute(risk = "Cold", dM, dF = dCp)) %>% mutate(risk = factor(risk, c("Heat", "Cold")))
lim <- c(-1.6, 3.2)
p3b <- ggplot(pl, aes(dM, pmin(pmax(dF, lim[1]), lim[2]), colour = risk)) + geom_abline(slope = -1, intercept = 0, linewidth = 0.35, linetype = "22", colour = "grey35") + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(size = 0.9, alpha = 0.7, stroke = 0) + facet_wrap(~ risk) + scale_colour_manual(values = risk_cols, guide = "none") + coord_cartesian(ylim = lim) +
  labs(x = "Change in log age standardised all cause death rate, 1990 to 2023", y = "Change in log age standardised\nattributable fraction", tag = "B")
save_fig(p3a / p3b + plot_layout(heights = c(1.15, 1)), "fig3_drivers", 7.2, 7.4)

# ---- Figure 4: causes -------------------------------------------------------------------------------
gw <- ro("global_causes.csv"); sh <- ro("shapley_share_causes.csv")
short <- c("Cardiovascular diseases" = "Cardiovascular", "Chronic respiratory diseases" = "Chronic respiratory", "Diabetes and kidney diseases" = "Diabetes and kidney", "Respiratory infections and tuberculosis" = "Respiratory infections",
           "Self-harm and interpersonal violence" = "Self-harm and violence", "Transport injuries" = "Transport injuries", "Unintentional injuries excluding forces of nature" = "Unintentional injuries")
co <- rev(unname(short))
kit <- function(p0, p1, w0, w1) tibble(within = 100 * (p1 - p0) * (w0 + w1) / 2, mix = 100 * (w1 - w0) * (p0 + p1) / 2)
k4 <- bind_rows(bind_cols(cause = gw$cause, risk = "Heat", kit(gw$ph_1990, gw$ph_2023, gw$w_1990, gw$w_2023)), bind_cols(cause = gw$cause, risk = "Cold", kit(gw$pc_1990, gw$pc_2023, gw$w_1990, gw$w_2023))) %>%
  pivot_longer(c(within, mix), names_to = "part", values_to = "v") %>% mutate(cause = factor(short[cause], co), risk = factor(risk, c("Heat", "Cold")), part = factor(part, c("within", "mix"), c("Change in the cause specific attributable fraction", "Change in the cause's share of all deaths")))
p4a <- ggplot(k4, aes(v, cause, fill = risk, alpha = part)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) + facet_wrap(~ risk) + scale_fill_manual(values = risk_cols, guide = "none") +
  scale_alpha_manual(values = c(1, 0.4), name = NULL) + labs(x = "Contribution to the change in the attributable fraction of all deaths, 1990 to 2023 (percentage points)", y = NULL, tag = "A") + guides(alpha = guide_legend(override.aes = list(fill = "grey30")))
s4 <- sh %>% mutate(cause = factor(short[cause], co), side = factor(ifelse(side == "heat", "Heat", "Cold"), c("Heat", "Cold")))
p4b <- ggplot(s4, aes(points, cause, fill = side)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) + scale_fill_manual(values = risk_cols, name = "Attributable deaths from", breaks = c("Heat", "Cold")) +
  labs(x = "Shapley contribution to the change in the global heat share, 1990 to 2023 (percentage points)", y = NULL, tag = "B")
save_fig(p4a / p4b + plot_layout(heights = c(1, 1), guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0)), "fig4_causes", 7.2, 6.6)

# ---- Figure 5: temperature ---------------------------------------------------------------------------
ZT <- ro("climate_axis.csv"); tr <- ro("country_trends.csv")
zw <- ZT %>% select(location_id, period, z, temp) %>% pivot_wider(names_from = period, values_from = c(z, temp))
lg <- function(p) log(p / (100 - p)); brk <- c(0.01, 0.1, 1, 10, 50, 90, 99)
fitl <- ZT %>% group_by(period) %>% group_modify(~ { m <- lm(z ~ temp, .x); tibble(temp = range(.x$temp), z = predict(m, tibble(temp = range(.x$temp)))) }) %>% mutate(period = ifelse(period == "p1", "1990 to 1994", "2016 to 2020"))
p5a <- ggplot(zw) + geom_segment(aes(x = temp_p1, y = z_p1, xend = temp_p2, yend = z_p2), colour = "grey45", linewidth = 0.3, arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_line(data = fitl, aes(temp, z, linetype = period), colour = INK, linewidth = 0.5) + scale_linetype_manual(values = c("22", "solid"), name = "Least squares line") +
  scale_y_continuous(breaks = lg(brk), labels = brk, limits = c(lg(0.001), lg(99.97))) + labs(x = "Population weighted mean annual temperature (°C)", y = "Heat share (%, logit scale)", tag = "A")
lim5 <- c(-1.5, 4)
p5b <- ggplot(tr, aes(pmin(pmax(imp_z, lim5[1]), lim5[2]), pmin(pmax(dz, lim5[1]), lim5[2]))) + geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = "22", colour = "grey35") + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(size = 0.9, alpha = 0.7, colour = "grey25", stroke = 0) + coord_cartesian(xlim = lim5, ylim = lim5) +
  labs(x = "Change implied by the country's own annual\ntemperature response and warming trend", y = "Observed trend in log heat to cold ratio,\n1990 to 2020", tag = "B")
save_fig(p5a + p5b + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "fig5_temperature", 7.2, 3.9)

# ---- Figure 6: shock years ---------------------------------------------------------------------------
gd <- ro("global_disaster_component.csv")
lab <- gd %>% arrange(-h729) %>% head(5) %>% mutate(label = sprintf("%d", year))
p6a <- ggplot(gd, aes(year)) + geom_col(aes(y = h729 / 1e3), fill = HEAT, width = 0.75) + geom_text(data = lab, aes(y = h729 / 1e3, label = label), vjust = -0.5, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = NULL, y = "Heat attributable deaths from exposure\nto forces of nature, global (thousands)", tag = "A")
cv <- bind_rows(lapply(list.files("data/raw", "multi_location_trend__r33[78]__c491__deaths__all__1990to2023", full.names = TRUE), function(f) read_csv(f, show_col_types = FALSE) %>% transmute(location_id, location_name, year, risk = ifelse(rei_id == 337, "Heat", "Cold"), paf = percent_combined)))
cases <- ro("shock_cases.csv") %>% filter(share_729 >= 15) %>% transmute(location_name, ev_year = year, event = sprintf("%s %d", sub(" \\(.*", "", location_name), year))
es <- cv %>% inner_join(cases, by = "location_name", relationship = "many-to-many") %>% mutate(k = year - ev_year) %>% filter(abs(k) <= 3) %>% mutate(risk = factor(risk, c("Heat", "Cold")))
p6b <- ggplot(es, aes(k, paf, group = event, colour = risk)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_line(linewidth = 0.4, alpha = 0.85) + geom_point(size = 0.8) +
  geom_text_repel(data = es %>% filter(k == 0, risk == "Heat"), aes(label = event), size = 2.2, colour = INK, nudge_x = 1.3, direction = "y", segment.size = 0.15, min.segment.length = 0) +
  facet_wrap(~ risk, scales = "free_y") + scale_colour_manual(values = risk_cols, guide = "none") + scale_x_continuous(breaks = -3:3) + labs(x = "Years from the shock year", y = "Cardiovascular deaths attributed\nto the risk (% of cardiovascular deaths)", tag = "B")
save_fig(p6a / p6b + plot_layout(heights = c(0.8, 1)) & theme(legend.position = "bottom"), "fig6_shock_years", 7.2, 6.2)

# ---- Supplementary figures -----------------------------------------------------------------------------
sig <- ro("sigma_series.csv"); ay <- ro("anomaly_years.csv")
s1a <- ggplot(sig, aes(year, sd_z)) + geom_line(linewidth = 0.55, colour = INK) + labs(x = NULL, y = "SD of log heat to cold ratio across countries", tag = "A") + scale_y_continuous(limits = c(2.4, 3.6))
s1b <- ggplot(sig, aes(year, med_z)) + geom_line(linewidth = 0.55, colour = INK) + labs(x = NULL, y = "Median log heat to cold ratio", tag = "B")
s1c <- ggplot(ay %>% filter(year <= 2019), aes(year, n_hot)) + geom_col(fill = HEAT, width = 0.75) + labs(x = NULL, y = "Countries with heat fraction more than\n2 robust SD above trend", tag = "C")
s1d <- ggplot(ay, aes(year, med_res)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 2019.5, linewidth = 0.3, linetype = "22", colour = "grey50") + geom_line(linewidth = 0.55, colour = HEAT) + labs(x = NULL, y = "Median residual of log heat fraction\nfrom the 1990 to 2019 trend", tag = "D")
save_fig((s1a + s1b) / (s1c + s1d), "figS1_dispersion_anomalies", 7.2, 5.8)
ga <- ro("global_age_pattern.csv") %>% mutate(age_group = factor(age_group, unique(age_group[order(match(age_group_id, c(2, 3, 388, 389, 238, 34, 6:20, 30, 31, 32, 235)))])))
gl <- ga %>% select(year, age_group, Heat = h_n, Cold = c_n) %>% pivot_longer(c(Heat, Cold), names_to = "risk", values_to = "n") %>% mutate(risk = factor(risk, c("Heat", "Cold")), year = factor(year))
s2 <- ggplot(gl, aes(age_group, n / 1e3, fill = risk, alpha = year)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) + geom_hline(yintercept = 0, linewidth = 0.3) + facet_wrap(~ risk, ncol = 1, scales = "free_y") + scale_fill_manual(values = risk_cols, guide = "none") +
  scale_alpha_manual(values = c(0.45, 1), name = NULL) + labs(x = "Age group", y = "Attributable deaths, global (thousands)") + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6.5)) + guides(alpha = guide_legend(override.aes = list(fill = "grey30")))
save_fig(s2, "figS2_age_pattern", 7.2, 5.6)
s3 <- ggplot(A %>% filter(location_id != 1) %>% mutate(region = factor(region, reg_order)), aes(year, share_n)) + geom_line(linewidth = 0.5, colour = INK) + facet_wrap(~ region, ncol = 4, labeller = label_wrap_gen(24)) + labs(x = NULL, y = "Heat share of temperature attributable deaths (%)") + scale_x_continuous(breaks = c(1990, 2005, 2020))
save_fig(s3, "figS3_region_share", 7.2, 4.2)
mc <- ro("country_share_mc.csv") %>% filter(is.finite(share23), share23 >= 0, share23 <= 100) %>% mutate(region = factor(c("Central Europe, Eastern Europe, and Central Asia" = "Central and eastern Europe, central Asia", "High-income" = "High income", "Latin America and Caribbean" = "Latin America and Caribbean",
  "North Africa and Middle East" = "North Africa and Middle East", "South Asia" = "South Asia", "Southeast Asia, East Asia, and Oceania" = "Southeast and east Asia, Oceania", "Sub-Saharan Africa" = "Sub-Saharan Africa")[super_region_name], reg_order[-1])) %>%
  group_by(region) %>% mutate(rk = rank(share23, ties.method = "first")) %>% ungroup()
s4f <- ggplot(mc, aes(rk, share23)) + geom_hline(yintercept = 50, linewidth = 0.3, linetype = "22", colour = "grey40") + geom_linerange(aes(ymin = pmax(lo23, 0), ymax = pmin(hi23, 100)), colour = "grey55", linewidth = 0.3) + geom_point(size = 0.6, colour = INK) +
  geom_point(aes(y = pmin(pmax(share90, 0), 100)), shape = 1, size = 0.7, colour = "grey45") + facet_wrap(~ region, ncol = 4, scales = "free_x", labeller = label_wrap_gen(24)) + labs(x = "Countries ranked by heat share in 2023 within super region", y = "Heat share (%)") + theme(axis.text.x = element_blank())
save_fig(s4f, "figS4_country_intervals", 7.2, 4.6)
cat("figures written:", paste(list.files("figures", "pdf$"), collapse = ", "), "\n")
