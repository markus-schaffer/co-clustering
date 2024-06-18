# Title: Visualize SHM data
#
# Purpose: Visualize the SHM data of all 4798 selected buildings. 
# This script produces Figure 2.
#
#
# Data files:
#   - data/01_shm_data.fst
#
# Function files:
#   - functions/plot_style.R
#
# Author: M. Schaffer 
# Contact details: msch@build.aau.dk


# Import packages --------------------------------------------------------

library(fst)
library(data.table)
library(purrr)
library(ggplot2)
library(lubridate)

source("functions/plot_style.R")
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Import Data -------------------------------------------------------------

equi_dt <- read_fst("data/01_shm_data.fst", as.data.table = TRUE)
equi_dt[, day := as.Date.POSIXct(time_rounded, tz = "Europe/Copenhagen")]


# Plot --------------------------------------------------------------------

# Calculate daily energy use in kWh and columns needed for plotting
daily_dt <- equi_dt[, .(daily_energy = sum(energy_area) / 1000), by = list(heat_meter_id, day)]
daily_dt[, week := isoweek(day)]
daily_dt[, month := month(day, label = TRUE, abbr = TRUE, locale = "English_United States")]
daily_dt[, pseudo_week := as.factor(rleid(week)), by = heat_meter_id]
daily_dt[, plot_day := mean(day), by = pseudo_week]

# needed to get axis always in English
Sys.setlocale(category = "LC_TIME", locale = "English_United States.1252")
p <- ggplot(daily_dt, aes(x = plot_day, y = daily_energy, group = pseudo_week)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, fill = "#DDDDDD") +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(daily_dt[, min(day)], daily_dt[, max(day)])) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", name = NULL, guide = guide_axis(angle = 90), expand = c(0.004, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), name = bquote("energy use - kWh/" ~ m^2 / day)) +
  theme_nice()

ggsave(
  filename = "plots/Figure_2.pdf",
  plot = p,
  device = cairo_pdf,
  width = 140,
  height = 80,
  units = "mm"
)
