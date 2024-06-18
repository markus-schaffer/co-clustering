# Title: Multinomial logistic regression with group lasso penalty
#
# Purpose: Multinomial logistic regression with group lasso penalty of the BCs
# against the clusters and analyses the result
# 
# This script produces Figure 10, Figure 11, Figure B.18 and the data for Table 3
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
library(msgl)
library(caret, include.only = c("createFolds", "upSample"))
library(furrr)
library(arm, include.only = "rescale")
library(ConfusionTableR)
library(ggplot2)
library(ggnewscale)
library(yardstick, include.only = c("mcc_vec"))

source("functions/plot_style.R")

plots <- list()


# Load Data ---------------------------------------------------------------

data <- list(
  bbr = fread("data/02_characteristics_bbr.csv"),
  epc = fread("data/02_characteristics_epc.csv")
)
data[["combined"]] <- merge(data$bbr, data$epc, by = "heat_meter_id", sort = FALSE)

cluster_result <- readRDS("data/03_clustering_reordered.RDS")
cluster_result <- cluster_result$row_clust

# Prepare data for classification -----------------------------------------

# Define categorical variables and one hot encode them
walk(data, ~ .x[, grep("*code", names(.x)) := map(.SD, ~ as.factor(.x)), .SDcols = grep("*code", names(.x))])
predictors <- map(data, ~ model.matrix(~ 0 + ., data = .x[, -"heat_meter_id"], contrasts.arg = lapply(.x[, grep("*code", names(.x)), with = FALSE], contrasts, contrasts = FALSE)))

# Scale predictors
norm_predictors <- map(predictors, ~ apply(.x, 2, arm::rescale))

# Define groups (factorial variables) for grouped lasso
fn_group <- function(data) {
  group <- gsub("(code).*", "\\1", colnames(data)) |>
    rle()
  rep(1:length(group$lengths), group$lengths)
}
group_id <- map(norm_predictors, ~ fn_group(.x))


# Multinomial logistic regression with group lasso penalty ----------------

# Get inner and outer folds for nested CV
set.seed(123)
outer_folds <- createFolds(cluster_result, list = TRUE, k = 5)
set.seed(123)
inner_folds <- map(outer_folds, ~ createFolds(cluster_result[-.x], k = 10))

# Set alpha to grouped lasso (0)
alpha <- 0

fn_lasso_optim <- function(x, y, outer_fold, inner_folds, group_id, data_name, alpha = 0, lambda = exp(seq(log(1), log(1e-5), length.out = 100)), nr_workers = 3) {
  cl <- makeCluster(nr_workers)
  registerDoParallel(cl)
  cv_fit <- msgl::cv(
    x = x[-outer_fold, ],
    classes = y[-outer_fold],
    grouping = group_id,
    cv.indices = inner_folds,
    lambda = lambda,
    alpha = alpha,
    standardize = FALSE,
    use_parallel = TRUE
  )
  stopCluster(cl)

  # Evaluate CV
  cv_mcc_mean <- apply(cv_fit$classes, 2, function(x) mcc_vec(truth = cv_fit$classes.true, estimate = factor(x, levels = as.character(unique(sort(cv_fit$classes.true))))))
  cv_best_fit <- data.table("mcc" = cv_mcc_mean[which.max(cv_mcc_mean)], "lambda" = cv_fit$lambda[[1]][[which.max(cv_mcc_mean)]])

  # Use best model parameters based on CV and fit model again as  msgl::cv does not return models
  # Two lambdas needed to force that the supplied lambda is used
  test_fit <- msgl::fit(
    x = x[-outer_fold, ],
    classes = y[-outer_fold],
    grouping = group_id,
    alpha = alpha,
    lambda = rep(cv_best_fit[, unique(lambda)], 2),
    standardize = FALSE
  )

  # MCC on test data
  test_pred <- predict(test_fit, x[outer_fold, ])
  test_mcc <- mcc_vec(truth = y[outer_fold], estimate = factor(test_pred$classes[, 1], levels = as.character(unique(sort(y[outer_fold])))))

  # Confusion Matrix
  test_confuse <- multi_class_cm(train_labels = factor(test_pred$classes[, 1], levels = as.character(unique(sort(y[outer_fold])))), truth_labels = y[outer_fold])

  return(
    list(
      "dataset" = data_name,
      "features" = features(test_fit)[[1]], # Selected features
      "test_mcc" = test_mcc,
      "test_confuse" = test_confuse,
      "cv" = cv_mcc_mean
    )
  )
}

plan(list(tweak(multisession, workers = 3), tweak(multisession, workers = 5)))
lasso_optim_result <- future_imap(norm_predictors, function(pred, pred_name) {
  future_imap(outer_folds, function(outer_fold, outer_fold_name) {
    fn_lasso_optim(
      x = pred,
      y = cluster_result,
      outer_fold = outer_fold,
      inner_folds = inner_folds[[outer_fold_name]],
      group_id = group_id[[pred_name]],
      data_name = pred_name,
      alpha = alpha,
      lambda = exp(seq(log(0.1), log(1e-5), length.out = 1e3)),
    )
  })
})
plan(sequential)

# saveRDS(lasso_optim_result, "data/05_mlrgl.RDS")
# lasso_optim_result <- readRDS("data/05_mlrgl.RDS")


# Plot results ------------------------------------------------------------

## Selected variables -----------------------------------------------------
sel_var <- map(lasso_optim_result, function(dataset) {
  imap(dataset, function(folds, folds_name) {
    data.table("variable" = folds$features, "fold" = folds_name)
  }) |> rbindlist()
})
sel_var <- rbindlist(sel_var, idcol = "dataset")
sel_var[, variable := gsub("(code).*", "\\1", variable)]
sel_var <- sel_var[, list(variable = unique(variable)), by = c("dataset", "fold")]
sel_var[dataset == "bbr", dataset := "BBR"]
sel_var[dataset == "epc", dataset := "EPC"]
sel_var[dataset == "combined", dataset := "Combined"]
sel_var[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]

fold_colors <- c(Fold1 = "#77AADD", Fold2 = "#44BB99", Fold3 = "#AAAA00", Fold4 = "#EEDD88", Fold5 = "#EE8866")

plots[["selected_vars"]] <-
  ggplot(sel_var, aes(x = variable, y = fold, fill = fold)) +
  geom_tile(color = "black", size = 0.01) +
  coord_flip() +
  facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = fold_colors) +
  scale_y_discrete(expand = c(0, 0), name = NULL) +
  scale_x_discrete(expand = c(0, 0), name = NULL) +
  theme_nice() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "plots/Figure_B_18.pdf",
  plot = plots[["selected_vars"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88 / (4 / 6),
  units = "mm"
)

## MCC and no. of selected variables Table ------------------------------------
test_mcc <- map(lasso_optim_result, function(dataset) {
  imap(dataset, function(folds, folds_name) {
    data.table("test_mcc" = folds$test_mcc, "fold" = folds_name)
  }) |> rbindlist()
})

test_mcc <- rbindlist(test_mcc, idcol = "dataset")
test_mcc <- test_mcc[, list(
  "mean_mcc" = mean(test_mcc),
  "min_mcc" = min(test_mcc),
  "max_mcc" = max(test_mcc)
),
by = dataset
]
test_mcc[dataset == "bbr", dataset := "BBR"]
test_mcc[dataset == "epc", dataset := "EPC"]
test_mcc[dataset == "combined", dataset := "Combined"]
test_mcc[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]
test_mcc <- merge(test_mcc, sel_var[, .("nr_var" = .N), by = c("dataset", "fold")][, .("nr_var" = mean(nr_var)), by = dataset], by = "dataset")
test_mcc[, nr_var:= nr_var-1] # Remove Intercept from variable count
print(test_mcc)


## Confusion matrix ---------------------------------------------------------
conf_matrix <- map(lasso_optim_result, function(data) {
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
conf_matrix[round(nomalised) == 0, fill := NA]

# Pretty up dataset names
conf_matrix[dataset == "bbr", dataset := "BBR"]
conf_matrix[dataset == "epc", dataset := "EPC"]
conf_matrix[dataset == "combined", dataset := "Combined"]
conf_matrix[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]

tp_colors <- c("#F7F7F7", "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837")
fp_fn_colors <- c("#F7F7F7", "#E7D4E8", "#C2A5CF", "#9970AB", "#762A83")

plots[["conf_matrix"]] <- ggplot() +
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
  filename = "plots/Figure_10.pdf",
  plot = plots[["conf_matrix"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88,
  units = "mm"
)


## MCC as function of the no of variables -----------------------------------
# This uses the MCC of the training data for simplicity
cv_mcc <- map(lasso_optim_result, function(dataset) {
  imap(dataset, function(folds, folds_name) {
    data.table("train_mcc" = folds$cv$mcc, "lambda" = folds$cv$lambda, "fold" = folds_name, nr_vars = folds$cv$nr_vars)
  }) |> rbindlist()
})
cv_mcc <- rbindlist(cv_mcc, idcol = "dataset")

# Pretty up dataset names
cv_mcc[dataset == "bbr", dataset := "BBR"]
cv_mcc[dataset == "epc", dataset := "EPC"]
cv_mcc[dataset == "combined", dataset := "Combined"]
cv_mcc[, dataset := factor(dataset, levels = c("BBR", "EPC", "Combined"))]

# Remove intercept from variable count
cv_mcc[, nr_vars := nr_vars - 1]

plots[["mcc_no_var"]] <- ggplot() +
  geom_line(data = cv_mcc,aes(x = nr_vars, y = train_mcc, group = fold, colour = fold), linewidth = 0.1 ) +
  geom_point(data = test_mcc, aes(x = nr_var, y = mean_mcc, fill = "MCC on test data")) +
  facet_grid(rows = vars(dataset)) +
  scale_fill_discrete(
    name = NULL,
    guide = guide_legend(override.aes = list(
      linetype = 0,
      shape = 16
    ), order = 2)
  ) +
  scale_colour_manual(
    values = fold_colors,
    name = NULL,
    breaks = "Fold1",
    labels = "MCC on train data",
    guide = guide_legend(
      order = 1,
      override.aes = list(
        color = "black"
      )
    )
  ) +
  scale_x_continuous(limits = c(1, 28), breaks = c(1, seq(5, 30, 5))) +
  labs(y = "MCC", x = "number of BCs") +
  theme_nice()

ggsave(
  filename = "plots/Figure_11.pdf",
  plot = plots[["mcc_no_var"]],
  device = cairo_pdf,
  dpi = 1200,
  width = 88,
  height = 88,
  units = "mm"
)
