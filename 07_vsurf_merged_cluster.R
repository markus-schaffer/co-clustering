# Title: VSURF + optimised RF with merged clusters
#
# Purpose: Runs VSURF and analyse result before optimising a random forest
# with the selected variables based on the merged cluster results.
#
# This script produces Figure 14, Figure 15 and data for Table 5
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
library(furrr)
library(VSURF)
library(caret, include.only = c("createFolds"))
library(ggplot2)
library(ConfusionTableR)
library(patchwork)
library(mlr)
library(tuneRanger)
library(ggnewscale)
library(yardstick, include.only = c("mcc_vec"))
library(forcats)


source("functions/plot_style.R")

if (!dir.exists(file.path("data/07_vsurf"))) {
  dir.create(file.path("data/07_vsurf"))
}

plots <- list()


# Import Data -------------------------------------------------------------

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
cluster_result <- fct_collapse(cluster_result,"E12" = c("E1", "E2"),"E34" = c("E3", "E4"),"E56" = c("E5", "E6"))
walk(data, ~ .x[, cluster := cluster_result])


# VSURF --------------------------------------------------------------------
# Outer folds of CV
set.seed(123)
outer_folds <- createFolds(cluster_result, list = TRUE, k = 5)

# VSURF leads to problems if executed in a function
for (i in 1:length(data)) {
  for (j in 1:length(outer_folds)) {
    print(paste(i, j))

    # Define training and test data
    train_X <- data[[i]][-outer_folds[[j]], -"cluster"] |> as.data.frame()
    train_Y <- data[[i]][-outer_folds[[j]], cluster]
    test_X <- data[[i]][outer_folds[[j]], -"cluster"] |> as.data.frame()
    test_Y <- data[[i]][outer_folds[[j]], cluster]

    set.seed(123)
    vsurf <-tryCatch(
        {VSURF(x = train_X, y = train_Y, parallel = TRUE, ncores = 50, verbose = FALSE,
                   ntree.thres = 2000,  nfor.thres = 50,
                   ntree.interp = 2000, nfor.interp = 25, 
                   ntree.pred = 2000,  nfor.pred = 25)
        }, error=function(r) {})
    
    if(inherits(vsurf, "error")) next
    
    pred <- predict(vsurf, newdata = test_X, step = c("thres", "interp", "pred"))
    # Test and confusion matrix
    test_mcc <- map_dbl(pred, ~ mean(.x != test_Y))
    test_mcc <- map_dbl(pred, ~mcc_vec(truth = test_Y, estimate = .x))
    test_confuse <- multi_class_cm(train_labels = pred$pred, truth_labels = test_Y)
    test_confuse_interp <- multi_class_cm(train_labels = pred$interp, truth_labels = test_Y)
    test_confuse_thres <- multi_class_cm(train_labels = pred$thres, truth_labels = test_Y)

    result <- data.table(
      "dataset" = names(data)[[i]],
      "test_mcc" = list(test_mcc),
      "test_confuse" = list(test_confuse),
      "test_confuse_interp" = list(test_confuse_interp),
      "test_confuse_thres" = list(test_confuse_thres),
      "variable_pred" = list(names(data[[i]][, -"cluster"])[vsurf$varselect.pred]),
      "variable_interp" = list(names(data[[i]][, -"cluster"])[vsurf$varselect.interp]),
      "variable_thres" = list(names(data[[i]][, -"cluster"])[vsurf$varselect.thres]),
      "fold" = names(outer_folds)[[j]]
    )

    saveRDS(result, paste0("data/07_vsurf/", i, "_", j, ".RDS"))
  }
}

# Collect result and delete the temp folder
files <- fs::dir_ls(path = "data/07_vsurf")
result_vsurf <- map(files, ~ readRDS(.x)) |> rbindlist()
saveRDS(result_vsurf, "data/07_vsurf_red.RDS")
unlink("data/06_vsurf", recursive = TRUE)

result_vsurf <- readRDS("data/07_vsurf_red.RDS")



# Analyse results of VSURF ---------------------------------------------------

## MCC and number of variables for each step of VSURF
test_mcc <- result_vsurf[, .(mcc = unlist(test_mcc), vsurf_step = names(unlist(test_mcc))), by = c("dataset", "fold")]
test_mcc <- test_mcc[, .(mean_mcc = mean(mcc), min_mcc = min(mcc), max_mcc = max(mcc)), by = c("dataset", "vsurf_step")]

nr_sel_var <- result_vsurf[, .(pred = length(unlist(variable_pred)), interp = length(unlist(variable_interp)), thres = length(unlist(variable_thres))), by = c("dataset", "fold")]
nr_sel_var <- melt(nr_sel_var, id.vars = c("dataset", "fold"), variable.name = "vsurf_step")
nr_sel_var <- nr_sel_var[, .(mean_var = mean(value), min_var = min(value), max_var = max(value)), by = c("dataset", "vsurf_step")]

merge(test_mcc, nr_sel_var, by = c("dataset", "vsurf_step")) |>
  print()

## Selected variables of predictor step --------------------------------------
sel_var <- result_vsurf[, .(variable_pred = unlist(variable_pred)), by = c("dataset", "fold")]
sel_var[, dataset := factor(dataset, levels = c("bbr", "epc", "combined"))]

# Pretty up dataset names
sel_var[dataset == "bbr", dataset := "BBR"]
sel_var[dataset == "epc", dataset := "EPC"]
sel_var[dataset == "combined", dataset := "Combined"]
sel_var[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]

fold_colors <- c("#77AADD", "#44BB99", "#AAAA00", "#EEDD88", "#EE8866")

plots[["selected_variables"]] <- ggplot(sel_var, aes(x = variable_pred, y = fold, fill = fold)) +
  geom_tile(color = "black", size = 0.01) +
  coord_flip() +
  facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = fold_colors) +
  scale_y_discrete(expand = c(0, 0), name = NULL) +
  scale_x_discrete(expand = c(0, 0), name = NULL) +
  theme_nice() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank()
  )

ggsave(
  filename = "plots/Figure_14.pdf",
  plot = plots[["selected_variables"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88 * 4 / 6,
  units = "mm"
)



# optimized RF --------------------------------------------------------------

# Define custom MCC measure for tuneRangeras the MCC of mlr (used under the
# hood) does not support multiclass
my_mcc_sub <- function(task, model, pred, feats, extra.args) {
  na_rm <- TRUE
  case_weights <- NULL
  estimator <- yardstick:::finalize_estimator(pred$data$truth, metric_class = "mcc")
  mcc_impl <- function(truth, response, ..., case_weights = NULL) {
    data <- yardstick:::yardstick_table(pred$data$truth, pred$data$response, case_weights = case_weights)
    yardstick:::mcc_table_impl(data, estimator)
  }
  yardstick:::metric_vec_template(
    metric_impl = mcc_impl, truth = pred$data$truth,
    estimate = pred$data$response, na_rm = na_rm, estimator = estimator,
    case_weights = case_weights, cls = "factor"
  )
}

my_mcc <- makeMeasure(
  id = "my.mcc", minimize = FALSE,
  properties = c("classif.multi", "classif"), fun = my_mcc_sub, best = 1, worst = -1
)

# Select only variables based on VSURF step
var <- list(bbr = c("renovation_code", "representative_year"), epc = c("vent_nat_winter", "total_transmission", "heating_temp_diff"), combined = c("vent_nat_winter", "total_transmission", "representative_year", "heating_temp_diff"))
walk2(data, var, ~ .x[, names(.x)[!names(.x) %in% c(.y, "cluster")] := NULL])


## Tune RF using tuneRanger --------------------------------------------------

fn_ranger <- function(data, outer_fold, data_name, threads = 3) {
  # Optimize hyperparameter of RF
  data_task <- makeClassifTask(data = as.data.frame(data[-outer_fold]), target = "cluster")
  res <- tuneRanger(data_task,
    measure = list(my_mcc), num.trees = 2000,
    num.threads = threads, iters = 70, save.file.path = NULL,
    tune.parameters = c(
      "mtry", "min.node.size",
      "sample.fraction"
    ),
    parameters = list(replace = TRUE, respect.unordered.factors = "order", "importance" = "none"),
    show.info = FALSE
  )

  data_pred <- predict(res$model, newdata = as.data.frame(data[outer_fold]))

  # Test error and confusion matrix
  test_mcc <- mcc_vec(truth = data_pred$data$truth, estimate = data_pred$data$response)
  test_confuse <- multi_class_cm(train_labels = data_pred$data$response, truth_labels = data_pred$data$truth)


  list(
    "dataset" = data_name,
    "test_mcc" = test_mcc,
    "test_confuse" = test_confuse
  )
}

plan(list(tweak(multisession, workers = 3), tweak(multisession, workers = 5)))
result_ranger <- future_imap(data, function(single_data, data_name) {
  future_map(outer_folds, function(outer_fold) {
    fn_ranger(
      data = single_data,
      outer_fold = outer_fold,
      data_name = data_name,
      threads = 3
    )
  }, .options = furrr_options(seed = set.seed(123)))
}, .options = furrr_options(seed = set.seed(123)))


saveRDS(result_ranger, "data/07_ranger.RDS")
# result_ranger <- readRDS("data/07_ranger.RDS")

## Analyse results of tuneRanger --------------------------------------------

# Calculate MCC
ranger_mcc <- map(result_ranger, function(data) {
  map(data, function(fold) {
    data.table(mcc = fold$test_mcc, dataset = fold$dataset)
  }) |> rbindlist()
}) |> rbindlist()

ranger_mcc <- ranger_mcc[, .(mean_mcc = mean(mcc), min_mcc = min(mcc), max_mcc = max(mcc)), by = c("dataset")]
ranger_mcc[, dataset := factor(dataset, levels = c("bbr", "epc", "combined"))]
print(ranger_mcc)

ranger_mcc[, map(.SD, ~ round(.x, 3)), .SDcols = is.numeric]


# Confusion Matrix
conf_matrix <- map(result_ranger, function(data) {
  imap(data, function(fold, fold_name) {
    cbind(as.data.table(fold$test_confuse$cm_tbl), data.table(dataset = fold$dataset, fold = fold_name))
  }) |> rbindlist()
}) |> rbindlist()
conf_matrix <- merge(conf_matrix, conf_matrix[, list(clas_sum = sum(N)), by = c("Reference", "dataset", "fold")], by = c("Reference", "dataset", "fold"))

# Normalize Confusion matrix and handle empty cases - i.e. with zero count
conf_matrix[, nomalised := N / clas_sum]
conf_matrix <- conf_matrix[, list(nomalised = mean(nomalised) * 100), by = c("Reference", "Prediction", "dataset")]
setnames(conf_matrix, "Reference", "Actual")
conf_matrix[, fill := nomalised]
conf_matrix[Actual != Prediction & round(nomalised) == 0, fill := NA]

# Pretty up dataset names
conf_matrix[dataset == "bbr", dataset := "BBR"]
conf_matrix[dataset == "epc", dataset := "EPC"]
conf_matrix[dataset == "combined", dataset := "Combined"]
conf_matrix[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]

tp_colors <- c("#F7F7F7", "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837")
fp_fn_colors <- c("#F7F7F7", "#E7D4E8", "#C2A5CF", "#9970AB", "#762A83")


plots[["conf_matrix_tune"]] <- ggplot() +
  geom_tile(data = conf_matrix[Actual == Prediction], aes(x = Prediction, y = Actual, fill = fill)) +
  scale_fill_gradientn(colors = (tp_colors), na.value = "#FFFFFF", name = "Correct", limits = c(0, 100), guide = guide_colorbar(order = 1)) +
  new_scale_fill() +
  geom_tile(data = conf_matrix[Actual != Prediction], aes(x = Prediction, y = Actual, fill = fill)) +
  scale_fill_gradientn(colors = (fp_fn_colors), na.value = "#FFFFFF", name = "Incorrect", limits = c(0, 100), guide = guide_colorbar(order = 2)) +
  geom_text(data = conf_matrix, aes(x = Prediction, y = Actual, label = round(nomalised)), color = "black", size = 2.1) +
  facet_grid(rows = vars(dataset)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_nice() +
  theme(legend.position = "right")

ggsave(
  filename = "plots/Figure_15.pdf",
  plot = plots[["conf_matrix_tune"]],
  device = cairo_pdf,
  width = 88,
  height = 88,
  units = "mm"
)
