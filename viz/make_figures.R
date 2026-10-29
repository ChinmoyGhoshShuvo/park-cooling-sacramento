# make_figures.R  -- figures and combined results table. Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
set.seed(SEED)
FIG <- file.path(DIR_OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
rd <- function(f) fread(file.path(DIR_OUT, f))
ring_mid <- c("0-50 m" = 25, "50-100 m" = 75, "100-150 m" = 125, "150-200 m" = 175, "200-300 m" = 250, "300-400 m" = 350, "400-500 m" = 450)
th <- theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
COL <- c("Headline: vs 400-500 m ring" = "#0072B2", "Spec: vs park-poor controls" = "#D55E00")

## 1. Love plots
lp <- function(B, title) {
  B <- B[ring != "park interior"]; B[, ring := factor(ring, levels = names(ring_mid))]
  L <- melt(B, id.vars = c("ring", "covariate"), measure.vars = c("smd_before", "smd_after"), variable.name = "when", value.name = "smd")
  L[, when := factor(when, c("smd_before", "smd_after"), c("Before matching", "After matching"))]
  ggplot(L, aes(abs(smd), covariate, shape = when, colour = when)) + geom_vline(xintercept = 0.1, linetype = 2, colour = "grey40") +
    geom_point(size = 2) + facet_wrap(~ring, nrow = 2) + scale_colour_manual(values = c("grey55", "#0072B2")) + scale_shape_manual(values = c(1, 16)) +
    labs(x = "|Standardized mean difference|  (dashed line = 0.1 target)", y = NULL, colour = NULL, shape = NULL, title = title) + th + theme(legend.position = "bottom")
}
ggsave(file.path(FIG, "love_plot_headline.png"), lp(rd("balance_table_headline.csv"), "Covariate balance, headline design (each ring vs 400-500 m ring)"), width = 11, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIG, "love_plot_spec.png"), lp(rd("balance_table.csv")[run == "built_main_exog"], "Covariate balance, original specification (vs controls >800 m from any green space)"), width = 11, height = 6, dpi = 300, bg = "white")

## 2. Map of pixels trimmed for lack of common support (original specification)
tm <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE); pk <- st_read(file.path(DIR_PROC, "parks_main_urban.gpkg"), quiet = TRUE)
cbd <- st_read(file.path(DIR_PROC, "cbd.gpkg"), quiet = TRUE)
Tr <- rd("../data_proc/trimmed_pixels.csv")[run == "built_main_exog"]
Tr <- unique(Tr, by = "cell"); xy <- xyFromCell(tm, Tr$cell); Tr[, `:=`(x = xy[, 1], y = xy[, 2])]
pool <- readRDS(file.path(DIR_PROC, "pool_locations.rds"))
nT <- uniqueN(pool[!is.na(ring_main) & bldg_dens150 >= 0.02 | interior_main == TRUE]$cell)
mp <- ggplot() + geom_sf(data = ua, fill = "grey97", colour = "grey40", linewidth = 0.3) + geom_sf(data = pk, fill = "#009E73", colour = NA, alpha = 0.5) +
  geom_point(data = Tr, aes(x, y, colour = trim_reason), size = 0.35, alpha = 0.7) + geom_sf(data = cbd, shape = 4, size = 3, stroke = 1.2) +
  scale_colour_manual(values = c("propensity above control range" = "#D55E00", "propensity below control range" = "#7B3294"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + coord_sf(datum = NA) +
  labs(title = "Treated pixels dropped for lack of common support (original specification)",
       subtitle = sprintf("%s of %s treated locations trimmed. Green = urban parks; X = CBD (State Capitol). Headline design trimmed only 14.", comma(nrow(Tr)), comma(nT)), x = NULL, y = NULL) +
  th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "map_trimmed_common_support.png"), mp, width = 8, height = 8.5, dpi = 300, bg = "white")
logmsg("Trimmed locations (original spec): ", nrow(Tr), " | within 5 km of CBD: ", sum(sqrt((Tr$x - st_coordinates(cbd)[1])^2 + (Tr$y - st_coordinates(cbd)[2])^2) < 5000))

## 3. Distance-decay curve
H <- rd("diagnostic_vs_400_500m_ring.csv")[, .(ring, att, lo = ci_lo, hi = ci_hi, design = "Headline: vs 400-500 m ring")]
S <- rd("results_att_by_ring.csv")[run == "built_main_exog" & ring != "park interior", .(ring, att = att_matched_diff, lo = ci95_lo_diff, hi = ci95_hi_diff, design = "Spec: vs park-poor controls")]
D <- rbind(H, S); D[, x := ring_mid[ring]]
dd <- ggplot(D, aes(x, att, colour = design, fill = design)) + geom_hline(yintercept = 0, colour = "grey30") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) + scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 300, 400, 500)) +
  labs(x = "Distance from park edge (m, ring midpoints)", y = "Effect on summer LST (deg C)", colour = NULL, fill = NULL,
       title = "Park spillover cooling, Sacramento (20 Landsat summer scenes, 2016-2025)",
       subtitle = "Ribbons: 95% CI, 1 km block bootstrap (500 reps). The 0-50 m ring is reported separately.") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "distance_decay_curve.png"), dd, width = 8.5, height = 5.5, dpi = 300, bg = "white")

## 4. Forest plot of per-park ATTs (50-150 m band, parks with >= 20 locations)
W <- rd("per_park_att.csv")
F <- W[n_loc_b50_150 >= 20][order(att_b50_150)]
F[, lab := sprintf("%s (%.1f ha)", unit_name, area_ha)]; F[, lab := factor(lab, levels = lab)]
fp <- ggplot(F, aes(att_b50_150, lab)) + geom_vline(xintercept = 0, colour = "grey30") + geom_errorbarh(aes(xmin = ci_lo_b50_150, xmax = ci_hi_b50_150), height = 0, colour = "grey50") +
  geom_point(aes(size = n_loc_b50_150), colour = "#0072B2") + scale_size_area(max_size = 3, name = "Locations") +
  labs(x = "Effect on summer LST 50-150 m from park edge (deg C, 95% CI)", y = NULL, title = "Per-park cooling beyond 50 m (headline design)",
       subtitle = sprintf("%d parks with at least 20 thinned locations in the 50-150 m band", nrow(F))) + th + theme(axis.text.y = element_text(size = 6.5))
ggsave(file.path(FIG, "forest_plot_per_park.png"), fp, width = 8.5, height = max(6, 0.14 * nrow(F) + 1.5), dpi = 300, bg = "white", limitsize = FALSE)

## 5. Park area vs cooling distance
Q <- W[n_loc_b0_50 >= 5 & n_loc_b50_150 >= 5 & n_loc_b150_200 >= 5]
AB <- rd("cooling_beyond_50m_by_park_area.csv")
sc <- ggplot(Q, aes(area_ha, cooling_distance_m)) + geom_jitter(width = 0, height = 8, alpha = 0.6, colour = "#0072B2") +
  scale_x_log10(labels = label_number(accuracy = 0.1)) + scale_y_continuous(breaks = c(0, 50, 150, 200), labels = c("none", "50 m", "150 m", "200 m")) +
  labs(x = "Park area (ha, log scale)", y = "Cooling distance", title = "Park area vs estimated cooling distance",
       subtitle = sprintf("%d parks with >=5 locations in every band; cooling distance = last contiguous band with 95%% CI < 0", nrow(Q))) + th
ab <- ggplot(AB, aes(factor(cls, levels = cls), att)) + geom_hline(yintercept = 0, colour = "grey30") + geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), colour = "#0072B2") +
  labs(x = "Park area class (ha)", y = "Effect 50-150 m (deg C)", title = "Pooled cooling beyond 50 m by park size") + th
ggsave(file.path(FIG, "park_area_vs_cooling_distance.png"), sc / ab, width = 8, height = 8, dpi = 300, bg = "white")

## 6. Combined results table
R7 <- rd("results_att_by_ring.csv"); NC <- rd("sensitivity_negative_control.csv"); RB <- rd("sensitivity_rosenbaum.csv")
on <- rd("headline_onsite_vs_400_500m.csv")
head <- rbind(rd("diagnostic_vs_400_500m_ring.csv")[, .(design = "HEADLINE: ring vs matched 400-500 m ring pixels", run = "built_main_exog", ring, method, att, se, ci_lo, ci_hi, att_dr, share_retained, max_abs_smd, unbalanced)],
              on[, .(design = "HEADLINE: ring vs matched 400-500 m ring pixels", run = "built_main_exog", ring, method = "NN", att, se, ci_lo, ci_hi, att_dr = NA, share_retained = NA, max_abs_smd = NA, unbalanced = "see balance_table_headline.csv")])
spec <- R7[, .(design = "SPEC: vs controls >800 m from any green space/water", run, ring, method, att = att_matched_diff, se = se_matched_diff, ci_lo = ci95_lo_diff, ci_hi = ci95_hi_diff,
               att_dr = att_doubly_robust, share_retained = share_treated_retained, max_abs_smd = max_abs_smd_after, unbalanced = unbalanced_covariates)]
out <- rbind(head, spec, fill = TRUE)
out <- merge(out, NC[, .(ring, run = "built_main_exog", negative_control = verdict)], by = c("run", "ring"), all.x = TRUE)
out <- merge(out, RB[, .(ring, run = "built_main_exog", rosenbaum_gamma = gamma_critical)], by = c("run", "ring"), all.x = TRUE)
out[, ring := factor(ring, levels = c(names(ring_mid), "park interior"))]; setorder(out, design, run, ring)
fwrite(out, file.path(DIR_OUT, "results_table.csv"))
logmsg("Figures: ", paste(list.files(FIG), collapse = ", "))
logmsg("STAGE_10_DONE")
