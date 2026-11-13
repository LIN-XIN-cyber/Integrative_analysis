args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "analysis/run_hcc_rf_pipeline.R"
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

local_lib <- Sys.getenv("HCC_R_LIB", unset = file.path(project_dir, "R_library"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

metadata_path <- if (length(args) >= 1) args[1] else
  "revise metadata final liver cancer add metagenomics.xlsx"
abundance_path <- if (length(args) >= 2) args[2] else
  "all_six_genus_level_tables.tsv"
output_dir <- if (length(args) >= 3) args[3] else file.path(project_dir, "outputs", "hcc_genus_rf")
excluded_cohorts <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",", fixed = TRUE)[[1]] else character()

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables_csv"), recursive = TRUE, showWarnings = FALSE)

required <- c("readxl", "readr", "dplyr", "tidyr", "ggplot2", "Boruta",
              "ranger", "pROC", "patchwork", "scales")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Boruta)
  library(ranger)
  library(pROC)
  library(patchwork)
})

source(file.path(project_dir, "analysis", "plot_style.R"))

SEED <- 20261107L
TRAIN_PROP <- 0.80
PREVALENCE_CUTOFF <- 0.10
MEAN_ABUNDANCE_CUTOFF <- 1e-4
BORUTA_P_VALUE <- 0.01
BORUTA_MAX_RUNS <- 300L
STABILITY_REPEATS <- 20L
STABILITY_SUBSAMPLE_PROP <- 0.50

STABILITY_CUTOFF <- 0.30
RF_TREES <- 1500L
CV_FOLDS <- 5L

set.seed(SEED)

theme_pub <- function(base_size = 11) {
  theme_hcc_publication(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, margin = margin(b = 7)),
      plot.subtitle = element_text(color = "#444444", margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      strip.text = element_text(face = "bold", color = "#111111")
    )
}

save_plot_both <- function(plot, stem, width = 8, height = 6) {
  ggsave(file.path(output_dir, "figures", paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(output_dir, "figures", paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
}

write_csv_out <- function(x, name) {
  readr::write_excel_csv(x, file.path(output_dir, "tables_csv", paste0(name, ".csv")), na = "")
}

message("Reading metadata and abundance data")
metadata_all <- readxl::read_excel(metadata_path, sheet = 1) |> as.data.frame()
needed_cols <- c("RunID", "ProjectID", "Status")
if (!all(needed_cols %in% names(metadata_all))) {
  stop("Metadata is missing required columns: ", paste(setdiff(needed_cols, names(metadata_all)), collapse = ", "))
}

metadata <- metadata_all |>
  mutate(
    RunID = as.character(RunID),
    ProjectID = factor(ProjectID),
    Status = factor(Status, levels = c("HD", "HCC"))
  ) |>
  filter(!is.na(RunID), !is.na(ProjectID), Status %in% c("HD", "HCC"))

excluded_by_design <- metadata |> filter(as.character(ProjectID) %in% excluded_cohorts)
metadata <- metadata |>
  filter(!as.character(ProjectID) %in% excluded_cohorts) |>
  mutate(ProjectID = droplevels(ProjectID))

if (anyDuplicated(metadata$RunID)) stop("Duplicate RunID values in metadata")

tax_raw <- readr::read_tsv(abundance_path, show_col_types = FALSE, progress = FALSE) |> as.data.frame()
names(tax_raw)[1] <- "taxon"
if (anyDuplicated(tax_raw$taxon)) stop("Duplicate taxon names in abundance table")
rownames(tax_raw) <- tax_raw$taxon
tax_raw$taxon <- NULL
tax_raw[] <- lapply(tax_raw, as.numeric)

matched_samples <- intersect(metadata$RunID, colnames(tax_raw))
excluded_metadata <- metadata |> filter(!RunID %in% matched_samples)
excluded_abundance <- setdiff(colnames(tax_raw), metadata$RunID)
metadata <- metadata |> filter(RunID %in% matched_samples)
tax_raw <- as.matrix(tax_raw[, metadata$RunID, drop = FALSE])
storage.mode(tax_raw) <- "double"

if (any(!is.finite(tax_raw)) || any(tax_raw < 0)) stop("Abundance matrix contains non-finite or negative values")
if (any(colSums(tax_raw) <= 0)) stop("At least one matched sample has zero total abundance")

# Remove only clearly non-genus placeholder rows; do not remove named uncultured genera automatically.
placeholder <- is.na(rownames(tax_raw)) | trimws(rownames(tax_raw)) %in% c("", "g__", "unclassified", "g__unclassified")
tax_raw <- tax_raw[!placeholder, , drop = FALSE]
tax_raw <- tax_raw[rowSums(tax_raw != 0) > 0, , drop = FALSE]
tax_tss <- sweep(tax_raw, 2, colSums(tax_raw), "/")

# Joint stratification retains every ProjectID x Status stratum in both partitions.
metadata$stratum <- interaction(metadata$ProjectID, metadata$Status, drop = TRUE, sep = "__")
set.seed(SEED)
split_tbl <- metadata |>
  group_by(stratum) |>
  group_modify(function(.x, .y) {
    n <- nrow(.x)
    if (n < 2) stop("A ProjectID x Status stratum has fewer than two samples: ", as.character(.y$stratum))
    n_train <- max(1L, min(n - 1L, round(TRAIN_PROP * n)))
    ids <- sample(seq_len(n), size = n, replace = FALSE)
    .x$Partition <- "Validation"
    .x$Partition[ids[seq_len(n_train)]] <- "Training"
    .x
  }) |>
  ungroup() |>
  arrange(ProjectID, Status, RunID)

metadata <- as.data.frame(split_tbl)
rownames(metadata) <- metadata$RunID
train_ids <- metadata$RunID[metadata$Partition == "Training"]
test_ids <- metadata$RunID[metadata$Partition == "Validation"]

split_balance <- metadata |>
  count(Partition, ProjectID, Status, name = "N") |>
  group_by(Partition, ProjectID) |>
  mutate(Proportion = N / sum(N)) |>
  ungroup()

# Filters are learned on training samples only, then held fixed for validation.
train_abd_all <- tax_tss[, train_ids, drop = FALSE]
prevalence <- rowMeans(train_abd_all > 0)
mean_abundance <- rowMeans(train_abd_all)
keep_feature <- prevalence > PREVALENCE_CUTOFF & mean_abundance > MEAN_ABUNDANCE_CUTOFF
if (!any(keep_feature)) stop("No taxa remain after prevalence/mean-abundance filtering")
tax_filt <- tax_tss[keep_feature, , drop = FALSE]
tax_filt <- sweep(tax_filt, 2, colSums(tax_filt), "/")

filter_table <- data.frame(
  Feature = rownames(tax_tss),
  TrainingPrevalence = prevalence,
  TrainingMeanAbundance = mean_abundance,
  PassPrevalence = prevalence > PREVALENCE_CUTOFF,
  PassMeanAbundance = mean_abundance > MEAN_ABUNDANCE_CUTOFF,
  Retained = keep_feature,
  stringsAsFactors = FALSE
) |>
  arrange(desc(Retained), desc(TrainingMeanAbundance))

source_mmuphin_adjust_batch <- function() {
  mmu_dir <- file.path(project_dir, "vendor", "MMUPHin", "R")
  files <- c("helpers.R", "helpers_adjust_batch.R", "adjust_batch.R")
  if (!all(file.exists(file.path(mmu_dir, files)))) stop("Vendored MMUPHin source files are missing")
  source(file.path(mmu_dir, "helpers.R"), local = .GlobalEnv)
  source(file.path(mmu_dir, "helpers_adjust_batch.R"), local = .GlobalEnv)
  control_adjust_batch <<- list(
    zero_inflation = TRUE,
    pseudo_count = NULL,
    diagnostic_plot = NULL,
    conv = 1e-4,
    maxit = 1000,
    verbose = TRUE
  )
  source(file.path(mmu_dir, "adjust_batch.R"), local = .GlobalEnv)
}

message("Running MMUPHin batch adjustment without outcome covariates")
source_mmuphin_adjust_batch()
mmu_meta <- as.data.frame(metadata[match(colnames(tax_filt), metadata$RunID), c("ProjectID"), drop = FALSE])
rownames(mmu_meta) <- colnames(tax_filt)
mmu_meta$ProjectID <- droplevels(factor(mmu_meta$ProjectID))
mmu <- adjust_batch(
  feature_abd = tax_filt,
  batch = "ProjectID",
  covariates = NULL,
  data = mmu_meta,
  control = list(zero_inflation = TRUE, diagnostic_plot = NULL, verbose = TRUE)
)
tax_mmu <- mmu$feature_abd_adj
tax_mmu[tax_mmu < 0 & tax_mmu > -1e-12] <- 0
if (any(tax_mmu < 0) || any(!is.finite(tax_mmu))) stop("Invalid MMUPHin-adjusted abundance values")
tax_mmu <- sweep(tax_mmu, 2, colSums(tax_mmu), "/")

clr_with_training_pseudo <- function(abundance, training_ids) {
  positive_train <- abundance[, training_ids, drop = FALSE]
  positive_train <- positive_train[positive_train > 0]
  if (!length(positive_train)) stop("No positive training abundances for CLR")
  pseudo <- min(positive_train) / 2
  log_x <- log(abundance + pseudo)
  clr <- sweep(log_x, 2, colMeans(log_x), "-")
  list(data = t(clr), pseudo_count = pseudo)
}

clr_before <- clr_with_training_pseudo(tax_filt, train_ids)
clr_after <- clr_with_training_pseudo(tax_mmu, train_ids)

stratified_folds <- function(meta, k, seed) {
  set.seed(seed)
  fold <- integer(nrow(meta))
  strata <- interaction(meta$ProjectID, meta$Status, drop = TRUE)
  for (s in levels(strata)) {
    idx <- which(strata == s)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

half_sample_by_stratum <- function(meta, prop, seed) {
  set.seed(seed)
  strata <- interaction(meta$ProjectID, meta$Status, drop = TRUE)
  unlist(lapply(split(seq_len(nrow(meta)), strata), function(idx) {
    n_take <- max(1L, min(length(idx), floor(length(idx) * prop)))
    sample(idx, n_take, replace = FALSE)
  }), use.names = FALSE)
}

run_boruta_stability <- function(x_train, y_train, meta_train, pipeline_name) {
  message("Boruta full-training fit: ", pipeline_name)
  set.seed(SEED + match(pipeline_name, c("Before MMUPHin", "After MMUPHin")) * 1000L)
  full_fit <- Boruta(
    x = as.data.frame(x_train), y = y_train,
    maxRuns = BORUTA_MAX_RUNS, pValue = BORUTA_P_VALUE,
    doTrace = 0, holdHistory = TRUE, mcAdj = TRUE
  )
  full_stats <- attStats(full_fit)
  full_stats$Feature <- rownames(full_stats)
  rownames(full_stats) <- NULL

  selection_matrix <- matrix(FALSE, nrow = ncol(x_train), ncol = STABILITY_REPEATS,
                             dimnames = list(colnames(x_train), paste0("Repeat", seq_len(STABILITY_REPEATS))))
  for (b in seq_len(STABILITY_REPEATS)) {
    idx <- half_sample_by_stratum(meta_train, STABILITY_SUBSAMPLE_PROP,
                                  SEED + 10000L + b + match(pipeline_name, c("Before MMUPHin", "After MMUPHin")) * 1000L)
    set.seed(SEED + 20000L + b + match(pipeline_name, c("Before MMUPHin", "After MMUPHin")) * 1000L)
    fit_b <- try(Boruta(
      x = as.data.frame(x_train[idx, , drop = FALSE]), y = droplevels(y_train[idx]),
      maxRuns = 150L, pValue = BORUTA_P_VALUE,
      doTrace = 0, holdHistory = FALSE, mcAdj = TRUE
    ), silent = TRUE)
    if (!inherits(fit_b, "try-error")) {
      dec <- as.character(fit_b$finalDecision)
      names(dec) <- names(fit_b$finalDecision)
      selection_matrix[names(dec)[dec == "Confirmed"], b] <- TRUE
    }
    message("  Stability repeat ", b, "/", STABILITY_REPEATS)
  }

  stability <- rowMeans(selection_matrix)
  full_stats$StabilityFrequency <- stability[full_stats$Feature]
  full_stats$PassStabilityCutoff <- full_stats$StabilityFrequency >= STABILITY_CUTOFF
  full_stats$Selected <- full_stats$decision == "Confirmed" & full_stats$PassStabilityCutoff
  full_stats$Pipeline <- pipeline_name
  full_stats <- full_stats |>
    select(Pipeline, Feature, decision, Selected, StabilityFrequency, PassStabilityCutoff,
           meanImp, medianImp, minImp, maxImp, normHits) |>
    arrange(desc(Selected), desc(StabilityFrequency), desc(medianImp))

  selected <- full_stats$Feature[full_stats$Selected]
  fallback <- FALSE
  if (length(selected) < 2) {
    selected <- full_stats$Feature[full_stats$decision == "Confirmed"]
    fallback <- TRUE
  }
  if (length(selected) < 2) stop("Fewer than two Boruta-confirmed features for ", pipeline_name)
  list(fit = full_fit, table = full_stats, selected = selected, fallback = fallback)
}

youden_threshold <- function(y, prob) {
  r <- pROC::roc(y, prob, levels = c("HD", "HCC"), direction = "<", quiet = TRUE)
  out <- pROC::coords(r, x = "best", best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  as.numeric(out[1, "threshold"])
}

classification_metrics <- function(y, prob, threshold) {
  truth <- factor(y, levels = c("HD", "HCC"))
  pred <- factor(ifelse(prob >= threshold, "HCC", "HD"), levels = c("HD", "HCC"))
  tp <- sum(pred == "HCC" & truth == "HCC")
  tn <- sum(pred == "HD" & truth == "HD")
  fp <- sum(pred == "HCC" & truth == "HD")
  fn <- sum(pred == "HD" & truth == "HCC")
  roc_obj <- pROC::roc(truth, prob, levels = c("HD", "HCC"), direction = "<", quiet = TRUE)
  ci_auc <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  data.frame(
    N = length(truth), HCC = sum(truth == "HCC"), HD = sum(truth == "HD"),
    AUC = as.numeric(pROC::auc(roc_obj)), AUC_CI_Lower = ci_auc[1], AUC_CI_Upper = ci_auc[3],
    Threshold = threshold,
    Accuracy = (tp + tn) / length(truth),
    Sensitivity = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
    Specificity = ifelse(tn + fp > 0, tn / (tn + fp), NA_real_),
    BalancedAccuracy = mean(c(ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
                              ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)), na.rm = TRUE),
    TP = tp, TN = tn, FP = fp, FN = fn
  )
}

tune_and_fit_rf <- function(x_train, y_train, meta_train, selected, pipeline_name) {
  x <- x_train[, selected, drop = FALSE]
  p <- ncol(x)
  fold <- stratified_folds(meta_train, CV_FOLDS,
                           SEED + 30000L + match(pipeline_name, c("Before MMUPHin", "After MMUPHin")) * 1000L)
  mtry_vals <- unique(pmax(1L, pmin(p, round(c(sqrt(p), p / 4, p / 3, p / 2)))))
  grid <- expand.grid(
    mtry = mtry_vals,
    min.node.size = c(1L, 5L, 10L),
    sample.fraction = c(0.70, 0.90),
    stringsAsFactors = FALSE
  )
  grid$MeanCV_AUC <- NA_real_
  grid$SD_CV_AUC <- NA_real_

  class_weights <- table(y_train)
  class_weights <- sum(class_weights) / (length(class_weights) * class_weights)

  for (g in seq_len(nrow(grid))) {
    fold_auc <- numeric(CV_FOLDS)
    for (f in seq_len(CV_FOLDS)) {
      idx_tr <- fold != f
      idx_va <- fold == f
      fit <- ranger(
        x = as.data.frame(x[idx_tr, , drop = FALSE]), y = droplevels(y_train[idx_tr]),
        probability = TRUE, num.trees = 700,
        mtry = grid$mtry[g], min.node.size = grid$min.node.size[g],
        sample.fraction = grid$sample.fraction[g], replace = FALSE,
        class.weights = class_weights, importance = "permutation",
        seed = SEED + g * 100L + f, num.threads = 1
      )
      pr <- predict(fit, data = as.data.frame(x[idx_va, , drop = FALSE]))$predictions[, "HCC"]
      fold_auc[f] <- as.numeric(pROC::auc(pROC::roc(y_train[idx_va], pr,
                                                    levels = c("HD", "HCC"), direction = "<", quiet = TRUE)))
    }
    grid$MeanCV_AUC[g] <- mean(fold_auc)
    grid$SD_CV_AUC[g] <- sd(fold_auc)
  }
  best <- grid |>
    arrange(desc(MeanCV_AUC), SD_CV_AUC, mtry, min.node.size, sample.fraction) |>
    slice(1)

  # Locked out-of-fold predictions define the operating threshold using training data only.
  oof <- rep(NA_real_, nrow(x))
  for (f in seq_len(CV_FOLDS)) {
    idx_tr <- fold != f
    idx_va <- fold == f
    fit <- ranger(
      x = as.data.frame(x[idx_tr, , drop = FALSE]), y = droplevels(y_train[idx_tr]),
      probability = TRUE, num.trees = RF_TREES,
      mtry = best$mtry, min.node.size = best$min.node.size,
      sample.fraction = best$sample.fraction, replace = FALSE,
      class.weights = class_weights, importance = "permutation",
      seed = SEED + 40000L + f, num.threads = 1
    )
    oof[idx_va] <- predict(fit, data = as.data.frame(x[idx_va, , drop = FALSE]))$predictions[, "HCC"]
  }
  threshold <- youden_threshold(y_train, oof)

  final_fit <- ranger(
    x = as.data.frame(x), y = y_train,
    probability = TRUE, num.trees = RF_TREES,
    mtry = best$mtry, min.node.size = best$min.node.size,
    sample.fraction = best$sample.fraction, replace = FALSE,
    class.weights = class_weights, importance = "permutation",
    seed = SEED + 50000L + match(pipeline_name, c("Before MMUPHin", "After MMUPHin")) * 1000L,
    num.threads = 1
  )

  list(fit = final_fit, grid = grid |> mutate(Pipeline = pipeline_name), best = best,
       oof = oof, threshold = threshold, fold = fold)
}

run_pipeline <- function(clr_obj, pipeline_name) {
  x_all <- clr_obj$data
  x_train <- x_all[train_ids, , drop = FALSE]
  x_test <- x_all[test_ids, , drop = FALSE]
  meta_train <- metadata[train_ids, , drop = FALSE]
  meta_test <- metadata[test_ids, , drop = FALSE]
  y_train <- droplevels(meta_train$Status)
  y_test <- droplevels(meta_test$Status)

  bor <- run_boruta_stability(x_train, y_train, meta_train, pipeline_name)
  message("Selected features for ", pipeline_name, ": ", length(bor$selected))
  rf <- tune_and_fit_rf(x_train, y_train, meta_train, bor$selected, pipeline_name)
  test_prob <- predict(rf$fit, data = as.data.frame(x_test[, bor$selected, drop = FALSE]))$predictions[, "HCC"]

  train_metrics <- classification_metrics(y_train, rf$oof, rf$threshold) |>
    mutate(Pipeline = pipeline_name, Dataset = "Training OOF")
  test_metrics <- classification_metrics(y_test, test_prob, rf$threshold) |>
    mutate(Pipeline = pipeline_name, Dataset = "Validation")

  pred <- data.frame(
    RunID = test_ids,
    ProjectID = as.character(meta_test$ProjectID),
    Status = as.character(y_test),
    Probability_HCC = test_prob,
    PredictedClass = ifelse(test_prob >= rf$threshold, "HCC", "HD"),
    Correct = ifelse(test_prob >= rf$threshold, "HCC", "HD") == as.character(y_test),
    Pipeline = pipeline_name,
    stringsAsFactors = FALSE
  )

  cohort_auc <- pred |>
    group_by(Pipeline, ProjectID) |>
    group_modify(function(.x, .y) {
      if (n_distinct(.x$Status) < 2) return(data.frame(N = nrow(.x), AUC = NA_real_, Lower = NA_real_, Upper = NA_real_))
      r <- pROC::roc(factor(.x$Status, levels = c("HD", "HCC")), .x$Probability_HCC,
                     levels = c("HD", "HCC"), direction = "<", quiet = TRUE)
      ci <- try(as.numeric(pROC::ci.auc(r, method = "delong")), silent = TRUE)
      if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_, NA_real_)
      data.frame(N = nrow(.x), AUC = as.numeric(pROC::auc(r)), Lower = ci[1], Upper = ci[3])
    }) |>
    ungroup()
  overall <- test_metrics |>
    transmute(Pipeline, ProjectID = "Overall", N, AUC, Lower = AUC_CI_Lower, Upper = AUC_CI_Upper)
  cohort_auc <- bind_rows(cohort_auc, overall)

  importance <- data.frame(
    Pipeline = pipeline_name,
    Feature = names(rf$fit$variable.importance),
    RF_PermutationImportance = as.numeric(rf$fit$variable.importance),
    stringsAsFactors = FALSE
  ) |>
    left_join(bor$table |> select(Feature, StabilityFrequency, medianImp), by = "Feature") |>
    arrange(desc(RF_PermutationImportance))

  list(
    pipeline = pipeline_name, pseudo_count = clr_obj$pseudo_count,
    boruta = bor, rf = rf, metrics = bind_rows(train_metrics, test_metrics),
    predictions = pred, cohort_auc = cohort_auc, importance = importance,
    roc = pROC::roc(y_test, test_prob, levels = c("HD", "HCC"), direction = "<", quiet = TRUE)
  )
}

result_before <- run_pipeline(clr_before, "Before MMUPHin")
result_after <- run_pipeline(clr_after, "After MMUPHin")
results <- list(result_before, result_after)

metrics <- bind_rows(lapply(results, `[[`, "metrics")) |>
  select(Pipeline, Dataset, everything())
predictions <- bind_rows(lapply(results, `[[`, "predictions"))
boruta_table <- bind_rows(lapply(results, function(z) z$boruta$table))
importance_table <- bind_rows(lapply(results, `[[`, "importance"))
cohort_auc <- bind_rows(lapply(results, `[[`, "cohort_auc"))
tuning_grid <- bind_rows(lapply(results, function(z) z$rf$grid)) |>
  select(Pipeline, everything())

best_params <- bind_rows(lapply(results, function(z) {
  data.frame(
    Pipeline = z$pipeline,
    SelectedFeatures = length(z$boruta$selected),
    StabilityFallbackUsed = z$boruta$fallback,
    CLR_PseudoCount = z$pseudo_count,
    mtry = z$rf$best$mtry,
    min.node.size = z$rf$best$min.node.size,
    sample.fraction = z$rf$best$sample.fraction,
    MeanCV_AUC = z$rf$best$MeanCV_AUC,
    SD_CV_AUC = z$rf$best$SD_CV_AUC,
    OOF_YoudenThreshold = z$rf$threshold
  )
}))

selected_features <- bind_rows(lapply(results, function(z) {
  z$importance |>
    filter(Feature %in% z$boruta$selected) |>
    mutate(SelectedRank = rank(-RF_PermutationImportance, ties.method = "first")) |>
    arrange(SelectedRank)
}))

split_table_out <- metadata |>
  select(RunID, ProjectID, Status, Partition, everything()) |>
  select(-stratum)

excluded_out <- bind_rows(
  data.frame(Source = rep("Metadata", nrow(excluded_by_design)), SampleID = excluded_by_design$RunID,
             Reason = rep(paste0("Excluded cohort sensitivity analysis: ", paste(excluded_cohorts, collapse = ",")), nrow(excluded_by_design)), stringsAsFactors = FALSE),
  data.frame(Source = rep("Metadata", nrow(excluded_metadata)), SampleID = excluded_metadata$RunID,
             Reason = rep("No matching abundance column", nrow(excluded_metadata)), stringsAsFactors = FALSE),
  data.frame(Source = rep("Abundance", length(excluded_abundance)), SampleID = excluded_abundance,
             Reason = rep("No matching metadata RunID", length(excluded_abundance)), stringsAsFactors = FALSE)
)

method_parameters <- data.frame(
  Parameter = c(
    "Random seed", "Training fraction", "Split stratification", "Prevalence cutoff",
    "Mean abundance cutoff", "Feature-filter fitting data", "Compositional transform",
    "CLR pseudocount", "MMUPHin batch", "MMUPHin covariates", "MMUPHin zero inflation",
    "Boruta p-value", "Boruta maximum runs", "Stability repeats", "Stability subsample fraction",
    "Candidate recurrence-frequency cutoff", "Random forest trees", "RF tuning", "Operating threshold",
    "Positive class", "Excluded cohorts"
  ),
  Value = c(
    SEED, TRAIN_PROP, "ProjectID x Status", PREVALENCE_CUTOFF,
    MEAN_ABUNDANCE_CUTOFF, "Training only", "TSS then CLR",
    "Half the minimum positive training abundance, separately per pipeline",
    "ProjectID", "None (Status not used)", TRUE,
    BORUTA_P_VALUE, BORUTA_MAX_RUNS, STABILITY_REPEATS, STABILITY_SUBSAMPLE_PROP,
    STABILITY_CUTOFF, RF_TREES, "5-fold ProjectID x Status stratified CV on training data",
    "Youden J from locked out-of-fold training predictions", "HCC",
    if (length(excluded_cohorts)) paste(excluded_cohorts, collapse = ",") else "None"
  ),
  stringsAsFactors = FALSE
)

session_lines <- capture.output(sessionInfo())
session_table <- data.frame(SessionInfo = session_lines, stringsAsFactors = FALSE)

write_csv_out(metrics, "model_metrics")
write_csv_out(best_params, "best_model_parameters")
write_csv_out(selected_features, "selected_features")
write_csv_out(boruta_table, "boruta_all_features")
write_csv_out(importance_table, "rf_feature_importance")
write_csv_out(predictions, "validation_predictions")
write_csv_out(cohort_auc, "validation_auc_by_cohort")
write_csv_out(tuning_grid, "rf_tuning_grid")
write_csv_out(split_table_out, "sample_split")
write_csv_out(split_balance, "split_balance")
write_csv_out(filter_table, "training_feature_filter")
write_csv_out(excluded_out, "excluded_samples")
write_csv_out(method_parameters, "method_parameters")
write_csv_out(session_table, "session_info")

# ROC plot with both locked holdout curves.
roc_df <- bind_rows(lapply(results, function(z) {
  data.frame(
    Specificity = z$roc$specificities,
    Sensitivity = z$roc$sensitivities,
    Pipeline = z$pipeline,
    stringsAsFactors = FALSE
  )
}))
auc_labels <- metrics |>
  filter(Dataset == "Validation") |>
  mutate(Label = sprintf("%s: AUC %.3f (95%% CI %.3f–%.3f)", Pipeline, AUC, AUC_CI_Lower, AUC_CI_Upper))

p_roc <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Pipeline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#9CA3AF", linewidth = 0.6) +
  geom_step(linewidth = 1.15) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = HCC_PIPELINE_COLORS) +
  labs(
    title = "HCC classification in the held-out 20% validation set",
    subtitle = paste(auc_labels$Label, collapse = "\n"),
    x = "1 − Specificity", y = "Sensitivity", color = NULL
  ) + theme_pub()
save_plot_both(p_roc, "roc_validation_comparison", 8.4, 6.6)

# Boruta and RF importance for selected features.
top_importance <- importance_table |>
  group_by(Pipeline) |>
  slice_max(RF_PermutationImportance, n = 20, with_ties = FALSE) |>
  ungroup() |>
  mutate(FeaturePlot = reorder(paste(Pipeline, Feature, sep = "___"), RF_PermutationImportance))
p_imp <- ggplot(top_importance, aes(x = RF_PermutationImportance, y = FeaturePlot, fill = Pipeline)) +
  geom_col(width = 0.72) +
  facet_wrap(~Pipeline, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_fill_manual(values = HCC_PIPELINE_COLORS) +
  labs(title = "Random forest importance of selected genera", x = "Permutation importance", y = NULL) +
  theme_pub() + theme(legend.position = "none", axis.text.y = element_text(size = 8.5))
save_plot_both(p_imp, "selected_feature_importance", 11.5, 7.5)

p_stability <- boruta_table |>
  filter(decision == "Confirmed" | StabilityFrequency >= 0.20) |>
  group_by(Pipeline) |>
  slice_max(StabilityFrequency, n = 25, with_ties = FALSE) |>
  ungroup() |>
  mutate(FeaturePlot = reorder(paste(Pipeline, Feature, sep = "___"), StabilityFrequency)) |>
  ggplot(aes(x = StabilityFrequency, y = FeaturePlot, color = Selected)) +
  geom_point(size = 2.4) +
  geom_vline(xintercept = STABILITY_CUTOFF, linetype = "dashed", color = "#6B7280") +
  facet_wrap(~Pipeline, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c(`TRUE` = "#047857", `FALSE` = "#9CA3AF")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Boruta feature-selection recurrence", subtitle = "Dashed line: 30% candidate recurrence-frequency cutoff",
       x = "Selection frequency across half-sample repeats", y = NULL, color = "Final selection") +
  theme_pub() + theme(axis.text.y = element_text(size = 8.5))
save_plot_both(p_stability, "boruta_stability", 11.5, 8.0)

# Cohort-level validation AUC forest plot.
cohort_auc_plot <- cohort_auc |>
  mutate(ProjectID = factor(ProjectID, levels = rev(c(sort(unique(as.character(metadata$ProjectID))), "Overall"))))
p_cohort <- ggplot(cohort_auc_plot, aes(x = AUC, y = ProjectID, color = Pipeline)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#9CA3AF") +
  geom_errorbarh(aes(xmin = pmax(0, Lower), xmax = pmin(1, Upper)), height = 0.16,
                 position = position_dodge(width = 0.45), linewidth = 0.7, na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.45), size = 2.8) +
  scale_color_manual(values = HCC_PIPELINE_COLORS) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Validation AUC by cohort", subtitle = "Error bars are DeLong 95% confidence intervals",
       x = "AUC", y = NULL, color = NULL) + theme_pub()
save_plot_both(p_cohort, "validation_auc_by_cohort", 8.8, 5.8)

# Split balance figure.
p_split <- split_balance |>
  ggplot(aes(x = ProjectID, y = N, fill = Status)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  facet_wrap(~Partition, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = HCC_GROUP_COLORS) +
  labs(title = "Project and outcome balance after the stratified 80:20 split",
       x = NULL, y = "Samples", fill = "Status") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_both(p_split, "split_balance", 9.0, 7.0)

# PCA shows cohort structure before and after correction using the same filtered genera.
pca_scores <- bind_rows(lapply(list("Before MMUPHin" = clr_before$data, "After MMUPHin" = clr_after$data), function(x) {
  pc <- prcomp(x, center = TRUE, scale. = FALSE)
  data.frame(RunID = rownames(x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
             Var1 = summary(pc)$importance[2, 1], Var2 = summary(pc)$importance[2, 2])
}), .id = "Pipeline") |>
  left_join(metadata |> select(RunID, ProjectID, Status, Partition), by = "RunID")

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = ProjectID, shape = Status)) +
  geom_point(alpha = 0.76, size = 1.8) +
  facet_wrap(~Pipeline, scales = "free", ncol = 2) +
  labs(title = "Cohort structure before and after MMUPHin correction",
       x = "Principal component 1", y = "Principal component 2", color = "Project", shape = "Status") +
  scale_color_manual(values = HCC_COHORT_COLORS, drop = FALSE) +
  theme_pub() + theme(legend.position = "right")
save_plot_both(p_pca, "pca_before_after_mmuphin", 11.0, 5.7)

# Compact summary for downstream workbook authoring.
summary_text <- c(
  paste0("Matched samples: ", nrow(metadata)),
  paste0("Training samples: ", length(train_ids)),
  paste0("Validation samples: ", length(test_ids)),
  paste0("Retained genera after training-only filter: ", sum(keep_feature)),
  paste0("Before MMUPHin selected genera: ", length(result_before$boruta$selected)),
  paste0("After MMUPHin selected genera: ", length(result_after$boruta$selected))
)
writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"), useBytes = TRUE)

saveRDS(list(
  metadata = metadata, filtered_abundance = tax_filt, mmuphin_abundance = tax_mmu,
  results = results, parameters = method_parameters
), file.path(output_dir, "analysis_objects.rds"), compress = "xz")

message("Analysis complete: ", normalizePath(output_dir, winslash = "/"))
