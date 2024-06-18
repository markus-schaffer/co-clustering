# Title: Analyse Building Characteristics in respect to clusters
#
# Purpose: Visually analyse the distribution of the Building Characteristics in
# respect to the clusters.
# This script produces Figure 9 and the supplementary material figures. 
#
#
# Data files:
#   - data/02_characteristics_bbr.csv
#   - data/02_characteristics_epc.csv
#   - data/03_clustering_reordered.RDS
#
# Function files:
#   - functions/plot_style.R
#
# Author: M. Schaffer Contact details: msch@build.aau.dk


# Import packages --------------------------------------------------------
library(data.table)
library(fst)
library(purrr)
library(ggridges)
library(ggplot2)
library(forcats)
library(patchwork)
library(gghighlight)

if (!dir.exists("plots/supplementary/bcs_vs_cluster")) {
  dir.create("plots/supplementary/bcs_vs_cluster", recursive = T)
}

source("functions/plot_style.R")

plots <- list()

# Import Data -------------------------------------------------------------

data <- list(
  bbr = fread("data/02_characteristics_bbr.csv"),
  epc = fread("data/02_characteristics_epc.csv")
)

cluster_result <- readRDS("data/03_clustering_reordered.RDS")

# Prepare data ------------------------------------------------------------

walk(data, ~ .x[, grep("*code", names(.x)) := map(.SD, ~ as.factor(.x)), .SDcols = grep("*code", names(.x))])
walk(data, ~ .x[, heat_meter_id := NULL])
walk(data, ~ .x[, cluster := cluster_result$row_clust])

# Calculate cluster normalized shares for variables with only a few distinct values
# Warning is because numeric variables such as "no_bathrooms" are converted to characters
bbr_discr <- suppressWarnings(
  melt(data$bbr[, grep("representative_year|developed_area_ratio", names(data[["bbr"]]), invert = TRUE), with = FALSE],
    id.vars = "cluster",
    value.factor = TRUE,
    variable.factor = FALSE
  )
)

# Reorder so that numeric values have "natural" order
fn_reorder <- function(level) {
  suppressWarnings(num <- as.numeric(level))
  c(sort(num), sort(level[is.na(num)]))
}
bbr_discr[, value := fct_relevel(value, fn_reorder)]
bbr_discr <- bbr_discr[, .N, by = list(variable, cluster, value)]
bbr_discr[, share := N / sum(N), by = list(cluster, variable)]

epc_discr <- melt(data$epc[, grep("*code|cluster", names(data[["epc"]]), invert = FALSE), with = FALSE],
  id.vars = "cluster",
  value.factor = TRUE,
  variable.factor = FALSE
)
epc_discr <- epc_discr[, .N, by = list(variable, cluster, value)]
epc_discr[, share := N / sum(N), by = list(cluster, variable)]


# Warning is because integer variables are converted to double
bbr_cont <- suppressWarnings(
  melt(data$bbr[, grep("representative_year|developed_area_ratio|cluster", names(data[["bbr"]])), with = FALSE],
    id.vars = "cluster",
    value.factor = FALSE,
    variable.factor = FALSE
  )
)

epc_cont <- melt(data$epc[, grep("*code", names(data[["epc"]]), invert = TRUE), with = FALSE],
  id.vars = "cluster",
  value.factor = FALSE,
  variable.factor = FALSE
)



# Publication plot -------------------------------------------------------

# Binwidth calculation: Freedman–Diaconis rule
fn_bw <- function(vec) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (all(is.wholenumber(vec))) {
    (2 * IQR(vec) / length(vec)^(1 / 3)) |> round()
  } else {
    (2 * IQR(vec) / length(vec)^(1 / 3)) |>
      formatC(format = "fg", digits = 2) |>
      as.numeric()
  }
}

cluster_colors <- c("#77AADD", "#FFAABB", "#BBCC33", "#EE8866", "#44BB99", "#EEDD88")

# after_stat(density) causes an error in gghighlight
plots[["rep_year"]] <- ggplot(data$bbr, aes(x = representative_year, fill = cluster)) +
  geom_histogram(aes(y = ..density..), binwidth = fn_bw(data$bbr$representative_year), boundary = 0) +
  gghighlight(unhighlighted_params = list(fill = "gray95"), label_key = cluster) +
  facet_grid(rows = vars(cluster)) +
  scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_y_continuous(breaks = c(0, 0.1, .2)) +
  scale_x_continuous(name = "representative year") +
  theme_nice()

plots[["total_trans"]] <- ggplot(data$epc, aes(x = total_transmission, fill = cluster)) +
  geom_histogram(aes(y = ..density..), binwidth = fn_bw(data$epc$total_transmission), boundary = 0) +
  gghighlight(unhighlighted_params = list(fill = "gray95"), label_key = cluster) +
  facet_grid(rows = vars(cluster)) +
  scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_y_continuous(breaks = c(0, 0.1, .2)) +
  scale_x_continuous(name = "total transmission losses") +
  theme_nice()

plots[["no_bathrooms"]] <- ggplot(bbr_discr[variable == "no_bathroom"], aes(x = value, y = share, fill = cluster)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_y_continuous() +
  labs(x = "number of\nbathrooms", y = "share per cluster") +
  scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
  theme_nice()

layout <- c(
  area(t = 1, l = 1, b = 6, r = 5),
  area(t = 1, l = 6, b = 6, r = 10),
  area(t = 1, l = 11, b = 3, r = 13),
  area(t = 4, l = 11, b = 6, r = 13)
)


plots[["publication"]] <- plots[["rep_year"]] + plots[["total_trans"]] + plots[["no_bathrooms"]] + guide_area() +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  guides(fill = guide_legend(ncol = 2)) &
  theme(legend.key.size = unit(5, "mm"))


ggsave(
  filename = "plots/Figure_9.pdf",
  plot = plots[["publication"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 160,
  height = 140 * 4 / 6,
  units = "mm"
)



# Supplementary plots -------------------------------------------------------

## BBR plots ---------------------------------------------------------------

# Discrete Characteristics
fn_discrete_plot <- function(data, sel_var) {
  ggplot(data[variable == sel_var], aes(x = value, y = share, fill = cluster)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_y_continuous() +
    labs(x = sel_var, y = "share per cluster") +
    scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
    theme_nice()
}
bbr_discr_plots <- map(sort(unique(bbr_discr$variable)), ~ fn_discrete_plot(data = bbr_discr, sel_var = .x))

plots[["bbr_disc_1"]] <- wrap_plots(bbr_discr_plots[1:4]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom")

plots[["bbr_disc_2"]] <- wrap_plots(bbr_discr_plots[5:8]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom")


ggsave(
  filename = "plots/supplementary/bcs_vs_cluster/bbr_1.pdf",
  plot = plots[["bbr_disc_1"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 160,
  height = 160,
  units = "mm"
)

ggsave(
  filename = "plots/supplementary/bcs_vs_cluster/bbr_2.pdf",
  plot = plots[["bbr_disc_2"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 160,
  height = 160,
  units = "mm"
)


# Continuous Characteristics - for some variables the IQR is 0. They used fixed 60 bins.
fn_cont_plot <- function(data, sel_var) {
  if (fn_bw(data[variable == sel_var, value]) == 0) {
    ggplot(data[variable == sel_var], aes(x = value, fill = cluster)) +
      geom_histogram(aes(y = ..density..), bins = 50, boundary = 0) +
      gghighlight(unhighlighted_params = list(fill = "gray95"), label_key = cluster) +
      facet_grid(rows = vars(cluster)) +
      scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
      scale_x_continuous(name = sel_var) +
      theme_nice()
  } else {
    ggplot(data[variable == sel_var], aes(x = value, fill = cluster)) +
      geom_histogram(aes(y = ..density..), binwidth = fn_bw(data[variable == sel_var, value]), boundary = 0) +
      gghighlight(unhighlighted_params = list(fill = "gray95"), label_key = cluster) +
      facet_grid(rows = vars(cluster)) +
      scale_fill_discrete(type = cluster_colors, name = "energy cluster") +
      scale_x_continuous(name = sel_var) +
      theme_nice()
  }
}

plots[["bbr_cont"]] <- map(sort(unique(bbr_cont$variable)), ~ fn_cont_plot(data = bbr_cont, sel_var = .x)) |>
  wrap_plots() +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom")


ggsave(
  filename = "plots/supplementary/bcs_vs_cluster/bbr_3.pdf",
  plot = plots[["bbr_cont"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 160,
  height = 160,
  units = "mm"
)

## EPC plots --------------------------------------------------------------

# Discrete
plots[["epc_disc"]] <- fn_discrete_plot(data = epc_discr, sel_var = "has_heat_pump_code")

# Continuous
epc_cont_plots <- map(sort(unique(epc_cont$variable)), ~ fn_cont_plot(data = epc_cont, sel_var = .x))
epc_cont_plots <- split(epc_cont_plots, ceiling(seq_along(epc_cont_plots) / 4))

# Add the one discrete plot to the last continuous one
epc_cont_plots[[4]] <- c(epc_cont_plots[[4]], list(plots[["epc_disc"]]))

epc_cont_plots <- map(epc_cont_plots, ~ wrap_plots(.x, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom"))

iwalk(epc_cont_plots, ~
  ggsave(
    filename = paste0("plots/supplementary/bcs_vs_cluster/epc_", .y, ".pdf"),
    plot = .x,
    device = cairo_pdf,
    dpi = 1200,
    width = 160,
    height = 160,
    units = "mm"
  ))
