# Title: Analyse the building characteristics (BCs)
#
# Purpose: Analyse the Pearson correlation coefficient of the BCs (Figure A_17)
# and calculate the Generalised Variance Inflation Factor (GVIF) (Figure 3)
#
#
# Data files:
#   - data/02_characteristics_bbr.csv
#   - data/02_characteristics_epc.csv
#
# Function files:
#   - functions/plot_style.R
#
# Author: M. Schaffer
# Contact details: msch@build.aau.dk


# Import packages --------------------------------------------------------
library(data.table)
library(fst)
library(purrr)
library(corrplot)
library(car)
library(ggplot2)

source("functions/plot_style.R")

# Load data ---------------------------------------------------------------

bbr_data <- fread("data/02_characteristics_bbr.csv")
epc_data <- fread("data/02_characteristics_epc.csv")


# Correlation plot --------------------------------------------------------

# Make all categorical variables to factors
bbr_data[, grep("*code", names(bbr_data)) := map(.SD, ~ as.factor(.x)), .SDcols = grep("*code", names(bbr_data))]
epc_data[, grep("*code", names(epc_data)) := map(.SD, ~ as.factor(.x)), .SDcols = grep("*code", names(epc_data))]

combined_data <- merge(bbr_data, epc_data, by = "heat_meter_id", sort = FALSE)
fn_colnames <- function(x) {
  colnames(x) <- paste0("_", colnames(x))
  x
}

# One hot encode data
predictors <- list(
  bbr = model.matrix(~ 0 + ., data = bbr_data[, -"heat_meter_id"], contrasts.arg = lapply(bbr_data[, grep("*code", names(bbr_data)), with = FALSE], function(x) fn_colnames(contrasts(x, contrasts = FALSE)))),
  epc = model.matrix(~ 0 + ., data = epc_data[, -"heat_meter_id"], contrasts.arg = lapply(epc_data[, grep("*code", names(epc_data)), with = FALSE], function(x) fn_colnames(contrasts(x, contrasts = FALSE)))),
  combined = model.matrix(~ 0 + ., data = combined_data[, -"heat_meter_id"], contrasts.arg = lapply(combined_data[, grep("*code", names(combined_data)), with = FALSE], function(x) fn_colnames(contrasts(x, contrasts = FALSE))))
)


## Correlation plots ----------------------------------------------------------
corr_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")

# Width and height are in inch
pdf(file = "plots/Figure_A_17_BBR.pdf", width = 88 * 3 / 25.4, height = 88 * 3 / 25.4)
corrplot(cor(predictors$bbr),
  method = "square",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.5 * 3,
  is.corr = TRUE,
  col = colorRampPalette(corr_colors)(20),
  order = "alphabet",
  cl.pos = "b",
  cl.cex = 0.7 * 3,
  cl.length = 5,
  addgrid.col = "gray60"
)
dev.off()

pdf(file = "plots/Figure_A_17_EPC.pdf", width = 88 * 3 / 25.4, height = 88 * 3 / 25.4)
corrplot(cor(predictors$epc),
  method = "square",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.5 * 3,
  is.corr = TRUE,
  col = colorRampPalette(corr_colors)(20),
  order = "alphabet",
  cl.pos = "b",
  cl.cex = 0.7 * 3,
  cl.length = 5,
  addgrid.col = "gray60"
)
dev.off()


pdf(file = "plots/Figure_A_17.pdf", width = 88 * 6 / 25.4, height = 88 * 6 / 25.4)
corrplot(cor(predictors$combined),
  method = "square",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.75 * 3,
  is.corr = TRUE,
  col = colorRampPalette(corr_colors)(20),
  order = "alphabet",
  cl.pos = "b",
  cl.cex = 0.7 * 3,
  cl.length = 5,
  addgrid.col = "gray60"
)
dev.off()


# GVIF --------------------------------------------------------------------

# Use of the heat_meter_id as dummy response needed because of the way the (G)VIF is implemented by car
lm_model <- lm(heat_meter_id ~ ., data = combined_data)
vif_lm <- as.data.table((vif(lm_model)[, 3, drop = FALSE])^2, keep.rownames = TRUE)
setnames(vif_lm, c("GVIF^(1/(2*Df))", "rn"), c("GVIF_lm", "predictor"))

p_gvif <- ggplot(vif_lm, aes(x = predictor, y = GVIF_lm)) +
  geom_col() +
  geom_hline(aes(yintercept = 5), linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = 10, y = 5, label = "used threshold", vjust = -0.5, size = 8 / .pt) +
  scale_y_continuous(
    name = expression(bgroup("(", GVIF^{
      (1 / (2 * "\u00D7" * DF))
    }, ")")^{
      2
    }),
    expand = c(0, 0), limits = c(0, 5.5)
  ) +
  scale_x_discrete(name = NULL, guide = guide_axis(angle = 90)) +
  theme_nice()

ggsave(
  filename = "plots/Figure_3.pdf",
  plot = p_gvif,
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 70,
  units = "mm"
)
