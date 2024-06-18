# Title: Decision tree to communicate cluster differentiating BCs
#
# Purpose: Decision trees to understand what BCs and how differentiate the found clusters.
#
# This script produces Figure 16 and Figures for the supplementary materials
#
#
# Data files:
#   - data/02_characteristics_bbr.csv
#   - data/02_characteristics_epc.csv
#   - data/03_clustering_reordered.RDS
#
# Function files:
#   - None
#
# Author: M. Schaffer Contact details: msch@build.aau.dk


# Import packages --------------------------------------------------------
library(data.table)
library(fst)
library(purrr)
library(rpart.plot)
library(forcats)
library(caret)
library(yardstick, include.only = c("mcc_vec"))
library(doParallel)

if (!dir.exists("plots/supplementary/cluster_visualisation")) {
  dir.create("plots/supplementary/cluster_visualisation", recursive = T)
}

# Load Cluster and BCs  ---------------------------------------------------

data <- list(
  bbr = fread("data/02_characteristics_bbr.csv"),
  epc = fread("data/02_characteristics_epc.csv")
)
data[["combined"]] <- merge(data$bbr, data$epc, by = "heat_meter_id", sort = FALSE)
cluster_result <- readRDS(paste0("data/03_clustering_reordered.RDS"))

# Prepare data ------------------------------------------------------------

walk(data, ~ .x[, grep("*code", names(.x)) := map(.SD, ~ as.factor(.x)), .SDcols = grep("*code", names(.x))])
walk(data, ~ .x[, heat_meter_id := NULL])

cluster_result <- cluster_result$row_clust
cluster_result <- fct_collapse(cluster_result, "E12" = c("E1", "E2"), "E34" = c("E3", "E4"), "E56" = c("E5", "E6"))
walk(data, ~ .x[, cluster := cluster_result])

var <- list(bbr = c("renovation_code", "representative_year"), epc = c("vent_nat_winter", "total_transmission", "heating_temp_diff"), combined = c("vent_nat_winter", "total_transmission", "representative_year", "heating_temp_diff"))
walk2(data, var, ~ .x[, names(.x)[!names(.x) %in% c(.y, "cluster")] := NULL])

# Abbreviate names that will be used - known from initial results
setnames(data[[1]], "representative_year", "year")
setnames(data[[1]], "renovation_code", "renovated")

# Decision trees ----------------------------------------------------------

## BBR -------------------------------------------------------------------

set.seed(123)
bbr_tree <- rpart(cluster ~ ., data = data[["bbr"]], method = "class", cp = 0, minsplit = 2, minbucket = 1)

printcp(bbr_tree)
plotcp(bbr_tree)

bbr_tree_pruned <- prune(bbr_tree, cp = 0.00137)
bbr_tree_pruned_pred <- predict(bbr_tree_pruned, data[["bbr"]], type = "class")
bbr_mcc <- mcc_vec(estimate = bbr_tree_pruned_pred, truth = data[[1]][, cluster])
bbr_acc <- mean(bbr_tree_pruned_pred == data[["bbr"]][, cluster])

print(paste("For BBR the MCC is:", round(bbr_mcc, 3), "and the accuracy is:", round(bbr_acc, 3)))

pdf(file = "plots/Figure_16.pdf", width = 120 / 25.4, height = 88 / 25.4)
rpart.plot(bbr_tree_pruned, yesno = 1, extra = 2, faclen = 0, varlen = 0, tweak = 1.2, roundint = TRUE, branch.lty = 2, type = 0, uniform = TRUE, legend.y = NA)
dev.off()

## EPC --------------------------------------------------------------------

set.seed(123)
epc_tree <- rpart(cluster ~ ., data = data[["epc"]], method = "class", cp = 0, minsplit = 2, minbucket = 1)
printcp(epc_tree)
plotcp(epc_tree)
epc_tree_pruned <- prune(epc_tree, cp = 0.00343456)
epc_tree_pruned_pred <- predict(epc_tree_pruned, data[["epc"]], type = "class")
epc_mcc <- mcc_vec(estimate = epc_tree_pruned_pred, truth = data[["epc"]][, cluster])
epc_acc <- mean(epc_tree_pruned_pred == data[["epc"]][, cluster])

pdf(file = "plots/supplementary/cluster_visualisation/epc.pdf", width = 88 / 25.4, height = 60 / 25.4)
rpart.plot(epc_tree_pruned,
  main = paste("For the EPC data, the MCC is:", round(epc_mcc, 3), "\nand the accuracy is:", round(epc_acc, 3)),
  yesno = 1, extra = 2, faclen = 0, varlen = 0, tweak = 1.2, roundint = TRUE, branch.lty = 2, type = 0, uniform = TRUE, legend.y = NA
)
dev.off()


## Combined  ---------------------------------------------------------------

set.seed(123)
combined_tree <- rpart(cluster ~ ., data = data[["combined"]][, c("cluster", "vent_nat_winter", "total_transmission", "representative_year")], method = "class", cp = 0, minsplit = 2, minbucket = 1)
printcp(combined_tree)
plotcp(combined_tree)

combined_tree_pruned <- prune(combined_tree, cp = 0.00325380)
combined_tree_pruned_pred <- predict(combined_tree_pruned, data[["combined"]][, c("cluster", "vent_nat_winter", "total_transmission", "representative_year")], type = "class")
combined_mcc <- mcc_vec(estimate = combined_tree_pruned_pred, truth = data[["combined"]][, cluster])
combined_acc <- mean(combined_tree_pruned_pred == data[["bbr"]][, cluster])


pdf(file = "plots/supplementary/cluster_visualisation/combined.pdf", width = 88 / 25.4, height = 60 / 25.4)
rpart.plot(combined_tree_pruned,
  main = paste("For the combined data, the MCC is:", round(combined_mcc, 3), "\nand the accuracy is:", round(combined_acc, 3)),
  yesno = 1, extra = 2, faclen = 0, varlen = 0, tweak = 1.2, roundint = TRUE, branch.lty = 2, type = 0, uniform = TRUE, legend.y = NA
)
dev.off()
