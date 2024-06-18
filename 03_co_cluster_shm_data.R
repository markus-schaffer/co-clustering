# Title: Co-cluster SHM data
#
# Purpose: Performs the co-clustering and the analysis of the co-clustering of
# the SHM data. First the optimal number of basis functions is determined. Then
# the data is clustered before the temporal and energy use clusters are
# analysed.
#
# This script produces Figure 4 to Figure 8.
#
#
# Data files:
#   - data/01_shm_data.fst
#   - data/00_weather_data.csv
#
# Function files:
#   - functions/plot_style.R
#   - functions/fn_fun_LBM.R
#   - functions/fn_funLBM_main.R
#
# Author: M. Schaffer Contact details: msch@build.aau.dk


# Import packages --------------------------------------------------------
library(fst)
library(data.table)
library(fda)
library(furrr)
library(purrr)
library(ggplot2)
library(plyr)
library(progressr)
library(funLBM)
library(scales)
library(mgsub)
library(lubridate)
library(patchwork)
library(grid)

options(future.globals.maxSize = 3500 * 1024^2) # Needed as the default size is too small for the SHM data

source("functions/plot_style.R")

# They are only used to be able to get more information back. They could be
# replaced to the standard functions again.
source("functions/fn_fun_LBM.R")
source("functions/fn_funLBM_main.R")

plots <- list()

if (!dir.exists("data/clustering_temp")) {
  dir.create("data/clustering_temp")
}


# Import Data -------------------------------------------------------------

shm_dt <- read_fst("data/01_shm_data.fst", as.data.table = TRUE)
setorder(shm_dt, heat_meter_id, time_rounded)

weather <- fread("data/00_weather_data.csv")

# Convert SHM data to matrix
timepoints <- 24
meter <- shm_dt[, uniqueN(heat_meter_id)]
days <- shm_dt[, .N] / timepoints / meter
energy_array <- array(shm_dt$energy_area, dim = c(timepoints, days, meter))
energy_array <- aperm(energy_array, c(3, 2, 1))
dimensions <- dim(energy_array)
print(dimensions)

# Find best number of basis functions -------------------------------------

# Evaluate basis function based GCV - one value per day
fourier_basis <- function(data, timepoints, nbasis) {
  basisobj <- create.fourier.basis(rangeval = range(1:timepoints), nbasis = nbasis)
  ys <- smooth.basis(argvals = 1:timepoints, y = data, fdParobj = basisobj)
  data.table(
    gcv = ys$gcv,
    nbasis = nbasis
  )
}

# Run evaluation in parallel to reduce time
plan(multisession, workers = 5)
result <- future_map(seq(3, (timepoints - 1), 2), function(n) {
  fourier_basis(data = t(apply(energy_array, 3, cbind)), timepoints = timepoints, nbasis = n)
}, .options = furrr_options(scheduling = 1L, seed = 123)) |> rbindlist()
plan(sequential)

# Mean result over all days
result_mean <- result[, .(mean_gcv = mean(gcv)), by = nbasis]

plots[["no_basis"]] <- ggplot(result_mean, aes(x = nbasis, y = mean_gcv)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  geom_point(data = result_mean[mean_gcv == min(mean_gcv)], size = 2) +
  scale_x_continuous(breaks = seq(3, (timepoints - 1), 4), minor_breaks = seq(5, (timepoints - 1), 4), name = "no. of basis functions") +
  scale_y_continuous(limits = c(20, 70), minor_breaks = NULL, name = "GCV") +
  coord_cartesian(xlim = c(3, 19)) +
  theme_nice()


ggsave(
  filename = "plots/Figure_4.pdf",
  plot = plots[["no_basis"]],
  device = cairo_pdf,
  width = 88,
  height = 88 * (4 / 6),
  units = "mm"
)



# Co-clustering -----------------------------------------------------------

## Run co-clustering ------------------------------------------------------

# All possible combinations of row (K) and column (L) cluster
param <- expand.grid(list(row_k = 2:12, col_l = 2:9))

# Error function needed as clustering fails for some combination of row and
# column cluster. clustering_temp is just a backup in case something goes wrong
fn_best_cluster <- function(input_data, K = 5, L = 5, nbasis = 7) {
  lbl <- tryCatch(
    {
      fun_LBM(input_data,
        K = K, L = L, nbasis = nbasis, maxit = 300, burn = 50,
        init = "funFEM", basis.name = "fourier", nbinit = 3, gibbs.it = 3
      )
    },
    error = function(e) {
      "error"
    }
  )
  saveRDS(data.table(lbl = list(lbl), K = K, L = L, nbasis = nbasis), paste0("data/clustering_temp/LMB_", K, "_", L, "_.RDS"))
  pg()
  data.table(lbl = list(lbl), K = K, L = L, nbasis = nbasis)
}

plan(multisession, workers = 18)
with_progress({
  pg <- progressor(steps = nrow(param))
  cluster_result <- future_map2(param$row_k, param$col_l, ~ fn_best_cluster(energy_array, K = .x, L = .y, nbasis = 7),
    .options = furrr_options(seed = 123, chunk_size = 1L)
  )
})
plan(sequential)
# saveRDS(cluster_result, "data/03_clustering_result.RDS")
# cluster_result <- readRDS("data/03_clustering_result.RDS")

# Visualise co-clustering -------------------------------------------------

cluster_result[map_lgl(lbl, ~ .x[[1]] != "error"), icl := map_dbl(lbl, ~ .x$icl)]
cluster_result[, L := factor(paste0("T", L), levels = paste0("T", sort(unique(L))))]
cluster_result[, K := factor(paste0("E", K), levels = paste0("E", sort(unique(K))))]

# Best result from clustering
best_cluster <- cluster_result[icl == max(icl, na.rm = TRUE), lbl][[1]]


## Evaluate the convergence of ICL ----------------------------------------
lbm_lst <- cluster_result[map_lgl(lbl, ~ .x[[1]] != "error"), lbl]
loglike_plot <- map(lbm_lst, ~ data.table("loglik" = .x$loglik, "K" = .x$K, "L" = .x$L, it = 1:length(.x$loglik))) |>
  rbindlist()
loglike_plot[, K_L := paste0(K, "_", L)]

plots[["icl_convergence"]] <- ggplot() +
  geom_line(data = loglike_plot[K != best_cluster$K & L != best_cluster$L], aes(x = it, y = loglik, group = K_L), linewidth = 0.2) +
  geom_line(data = loglike_plot[K == best_cluster$K & L == best_cluster$L], aes(x = it, y = loglik, group = K_L), linewidth = 0.5, colour = "#DD3D2D") +
  scale_y_continuous(labels = label_scientific(digits = 3)) +
  labs(x = "iterations", y = "Complete log−likelihood") +
  theme_nice()


## Maximum ICL plot -------------------------------------------------------

icl_colors <- c(
  "#DDECBF", "#D0E7CA", "#C2E3D2", "#B5DDD8", "#A8D8DC",
  "#9BD2E1", "#8DCBE4", "#81C4E7", "#7BBCE7", "#7EB2E4", "#88A5DD", "#9398D2",
  "#9B8AC4", "#9D7DB2"
) |>
  rev()

plots[["icl_selection"]] <- ggplot(cluster_result, aes(x = L, y = K, fill = icl)) +
  geom_tile() +
  geom_text(data = cluster_result[icl == max(icl, na.rm = TRUE)], aes(x = L, y = K, label = "max", color = "NA"), size = 6 / .pt) +
  scale_x_discrete(name = "time cluster", expand = c(0, 0)) +
  scale_y_discrete(name = "energy use cluster", expand = c(0, 0)) +
  scale_fill_gradientn(colors = icl_colors, na.value = "#FFFFFF", labels = scientific, name = "ICL") +
  scale_color_manual(values = "black", name = "Failed", labels = NULL) +
  guides(color = guide_legend(override.aes = list(fill = "white", color = "white"), order = 2)) +
  theme_nice() +
  theme(
    legend.position = "right",
    legend.key = element_rect(colour = "black"),
    legend.box.margin = margin(0, 0, -15, 0)
  )

ggsave(
  filename = "plots/Figure_5.pdf",
  plot = plots[["icl_selection"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88 * 4 / 6,
  units = "mm"
)


## Rename cluster for more "natural" naming based on results ---------------

# This requires an initial analysis of the results beforehand

# temporal clusters
tc_org <- c("T3", "T2", "T4", "T6", "T1", "T5")
tc_new <- c("T6", "T5", "T4", "T3", "T2", "T1")

# Energy clusters
ec_org <- c("E5", "E2", "E1", "E4", "E6", "E3")
ec_new <- c("E6", "E5", "E4", "E3", "E2", "E1")

best_cluster$col_clust <- paste0("T", best_cluster$col_clust)
best_cluster$col_clust <- mgsub(best_cluster$col_clust, pattern = tc_org, replacement = tc_new) |>
  factor()
best_cluster$row_clust <- paste0("E", best_cluster$row_clust)
best_cluster$row_clust <- mgsub(best_cluster$row_clust, pattern = ec_org, replacement = ec_new) |>
  factor()

saveRDS(best_cluster, "data/03_clustering_reordered.RDS")


## Plot temporal clusters ---------------------------------------------------

# Get the days corresponding to the different temporal clusters
shm_dt[, col_clust := rep(rep(best_cluster$col_clust, each = timepoints), meter)]
shm_dt[, row_clust := rep(best_cluster$row_clust, each = days * timepoints)]
date_cluster <- data.table(
  day = with_tz(shm_dt$time_rounded[1] + days(seq(0, days - 1)), "CET"),
  cluster = best_cluster$col_clust
)

# Add weather data
date_cluster <- merge(date_cluster, weather, by = "day")

# Get weeks, and weekdays for plotting
date_cluster[, week := isoweek(day)]
date_cluster[, weekday := wday(day, label = TRUE, week_start = 1, locale = "English_United States")]
date_cluster[, month := month(day, label = TRUE, abbr = TRUE, locale = "English_United States")]
date_cluster[month == "Jan" & week > 10, week := 0]
date_cluster[, year := year(day)]
date_cluster[, pseudo_week := rleid(week), by = list(year, month)]

season_colors <- c("#1965B0", "#7BAFDE", "#4EB265", "#CAE0AB", "#F1932D", "#DC050C")
temp_colors <- c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")

plots[["temporal_clusters"]] <- ggplot(date_cluster, aes(x = pseudo_week, y = weekday, fill = cluster)) +
  geom_tile() +
  facet_grid(rows = vars(year), cols = vars(month), scale = "free") +
  scale_fill_manual(values = season_colors) +
  theme() +
  scale_x_continuous(expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_discrete(expand = c(0, 0), name = NULL) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "mm"),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "white", linewidth = 0.1),
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 8, colour = "black"),
    line = element_line(colour = "black", linewidth = 0.1),
    rect = element_rect(colour = "black", linewidth = 0.1),
    axis.ticks = element_line(colour = "black", linewidth = 0.1),
    legend.position = "none"
  )


plots[["temporal_temperature"]] <- ggplot(date_cluster, aes(x = pseudo_week, y = weekday, fill = day_mean_temp)) +
  geom_tile() +
  facet_grid(rows = vars(year), cols = vars(month), scale = "free") +
  scale_fill_gradientn(
    colors = temp_colors,
    rescaler = ~ rescale_mid(.x, mid = 0),
    name = "daily mean\nexterior T - °C",
    breaks = seq(-10, 50, 5)
  ) +
  guides(fill = guide_colourbar(barwidth = 10)) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_discrete(expand = c(0, 0), name = NULL) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "mm"),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "white", linewidth = 0.1),
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 8, colour = "black"),
    line = element_line(colour = "black", linewidth = 0.1),
    rect = element_rect(colour = "black", linewidth = 0.1),
    axis.ticks = element_line(colour = "black", linewidth = 0.1),
    legend.position = "bottom"
  )


date_cluster[, temperature_median := median(day_mean_temp), by = cluster]

plots[["temperature_clusters"]] <- ggplot(date_cluster, aes(x = reorder(cluster, -temperature_median), y = day_mean_temp, fill = cluster)) +
  geom_violin(color = "black", size = 0.1) +
  stat_summary(
    fun.min = median, fun = median, fun.max = median,
    geom = "errorbar",
    width = 0.45,
    linewidth = 0.75,
    aes(colour = "Mean"),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = season_colors, name = "Time cluster") +
  scale_color_manual(values = "black", name = "") +
  scale_x_discrete(breaks = NULL, name = NULL) +
  scale_y_continuous(limits = c(-10, NA), name = "daily mean exterior T - °C") +
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  theme_nice()


layout <- c(
  area(t = 1, l = 1, b = 4, r = 4),
  area(t = 5, l = 1, b = 8, r = 4),
  area(t = 9, l = 1, b = 10, r = 2),
  area(t = 9, l = 3, b = 10, r = 4)
)


plots[["temporal_publication"]] <- plots[["temporal_clusters"]] + plots[["temporal_temperature"]] + plots[["temperature_clusters"]] + guide_area() +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave(
  filename = "plots/Figure_6.pdf",
  plot = plots[["temporal_publication"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 180,
  height = 200,
  units = "mm"
)


### Solar irradiation --------------------------------------------------------

# Pearson Correlation Coefficient of temperature and solar irradiation
date_cluster[, cor(day_mean_rad, day_mean_temp)] |> round(3)

plots[["temporal_irrad"]] <- ggplot(date_cluster, aes(x = pseudo_week, y = weekday, fill = day_mean_rad)) +
  geom_tile() +
  facet_grid(rows = vars(year), cols = vars(month), scale = "free") +
  scale_fill_gradientn(
    colors = temp_colors,
    name = "daily mean\nglobal radiation - W/m²"
  ) +
  guides(fill = guide_colourbar(barwidth = 10)) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL, name = NULL) +
  scale_y_discrete(expand = c(0, 0), name = NULL) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "mm"),
    panel.spacing.y = unit(1, "mm"),
    strip.background = element_rect(fill = "white", linewidth = 0.1),
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 8, colour = "black"),
    line = element_line(colour = "black", linewidth = 0.1),
    rect = element_rect(colour = "black", linewidth = 0.1),
    axis.ticks = element_line(colour = "black", linewidth = 0.1),
    legend.position = "bottom"
  )

date_cluster[, irrad_median := median(day_mean_rad), by = cluster]
plots[["irrad_clusters"]] <- ggplot(date_cluster, aes(x = reorder(cluster, -day_mean_rad), y = day_mean_rad, fill = cluster)) +
  geom_violin(color = "black", size = 0.1) +
  stat_summary(
    fun.min = median, fun = median, fun.max = median,
    geom = "errorbar",
    width = 0.45,
    linewidth = 0.75,
    aes(colour = "Mean"),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = season_colors, name = "Time cluster") +
  scale_color_manual(values = "black", name = "") +
  scale_x_discrete(breaks = NULL, name = NULL) +
  scale_y_continuous(limits = c(-10, NA), name = "daily mean exterior T - °C") +
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  theme_nice()

plots[["temporal_irrad_final"]] <- plots[["temporal_clusters"]] + plots[["temporal_irrad"]] + plots[["irrad_clusters"]] + guide_area() +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")")


ggsave(
  filename = "plots/Figure_6_not_shown.pdf",
  plot = plots[["temporal_irrad_final"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 180,
  height = 200,
  units = "mm"
)


## Plot cluster means -------------------------------------------------------

# Evaluate the cluster mean functions as estimated by FunLBM to get discrete
# plot-able values

# https://stackoverflow.com/questions/63745186/how-can-i-plot-an-fda-object-using-ggplot2
TT <- best_cluster$T
K <- best_cluster$K
L <- best_cluster$L
grid <- 24
basis <- create.fourier.basis(c(0, TT), best_cluster$nbasis)
obj <- list(basis = basis, coefs = c(), fdnames = list(time = 1:TT, reps = c(), values = c()))
class(obj) <- "fd"
cluster_mean <- data.table()
for (k in 1:K) {
  for (l in 1:L) {
    xx <- seq(0, TT, len = grid * TT)
    obj$coefs <- replicate(2, best_cluster$prms$mu[k, l, ])
    obj_mat <- eval.fd(xx, obj)
    obj_df <- obj_mat |>
      as.data.table()
    obj_df[, V2 := NULL]
    obj_df[, x := xx]
    obj_df[, row := k]
    obj_df[, col := l]
    cluster_mean <- rbind(cluster_mean, obj_df)
  }
}
setnames(cluster_mean, "V1", "basis_eval")

# Rename to match used naming
cluster_mean[, row := paste0("E", row)]
cluster_mean[, row := mgsub(row, pattern = ec_org, replacement = ec_new) |>
  factor()]
cluster_mean[, col := paste0("T", col)]
cluster_mean[, col := mgsub(col, pattern = tc_org, replacement = tc_new) |>
  factor()]


plots[["cluster_mean"]] <- ggplot(cluster_mean, aes(x = x, y = basis_eval)) +
  facet_grid(row ~ col) +
  geom_line() +
  labs(y = bquote("cluster mean - Wh/" ~ m^2), x = "hour of the day") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, length.out = 4)) +
  theme_nice()


#  Recolor the temporal cluster backgrounds to match used colors
# https://stackoverflow.com/a/53457599/16697029

# Generate the ggplot2 plot grob
g <- grid.force(ggplotGrob(plots[["cluster_mean"]]))
# Get the names of grobs and their gPaths into a data.frame structure
grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
# Build optimal gPaths that will be later used to identify grobs and edit them
grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
grobs_df$gPath_full <- gsub(
  pattern = "layout::",
  replacement = "",
  x = grobs_df$gPath_full,
  fixed = TRUE
)
# Get the gPaths of the strip background grobs
strip_bg_gpath <- grobs_df$gPath_full[grepl(
  pattern = ".*strip\\.background.*",
  x = grobs_df$gPath_full
)]
# Get the gPaths of the strip titles
strip_txt_gpath <- grobs_df$gPath_full[grepl(
  pattern = "strip.*titleGrob.*text.*",
  x = grobs_df$gPath_full
)]
# Define the colors as seasonal colors plus white for rermaining fills
n_cols <- length(strip_bg_gpath)
fills <- c(season_colors, rep("#FFFFFF", n_cols - length(season_colors)))

# Edit the grobs - recolor them
for (i in 1:length(strip_bg_gpath)) {
  g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
}
plots[["cluster_mean"]] <- g


ggsave(
  filename = "plots/Figure_7.pdf",
  plot = plots[["cluster_mean"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 140,
  height = 140,
  units = "mm"
)


## Plot cluster proportion -------------------------------------------------

# Get proportion from results
cluster_proportions <- data.table(
  parameter = c(
    paste0("E", seq_along(best_cluster$prms$alpha)),
    paste0("T", seq_along(best_cluster$prms$beta))
  ),
  estimate = c(best_cluster$prms$alpha, best_cluster$prms$beta),
  group = c(
    rep("energy use cluster", length(best_cluster$prms$alpha)),
    rep("time cluster", length(best_cluster$prms$beta))
  )
)

# Rename the cluster names to be consistent
cluster_proportions[, parameter := mgsub(parameter, pattern = tc_org, replacement = tc_new)]
cluster_proportions[, parameter := mgsub(parameter, pattern = ec_org, replacement = ec_new)]

plots[["cluster_share"]] <- ggplot(cluster_proportions, aes(x = parameter, y = estimate)) +
  geom_col(fill = "gray80", color = "black", size = 0.05) +
  facet_grid(cols = vars(group), scale = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.3)) +
  labs(x = "cluster", y = "proportion") +
  theme_nice()

ggsave(
  filename = "plots/Figure_8.pdf",
  plot = plots[["cluster_share"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88 * 4 / 6,
  units = "mm"
)
