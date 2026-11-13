project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
local_lib <- Sys.getenv("HCC_R_LIB", unset = file.path(project_dir, "R_library"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source(file.path(project_dir, "analysis", "plot_style.R"))

out_dir <- file.path(project_dir, "outputs", "reviewer_16s_only_sensitivity")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

shap6 <- readRDS(file.path(
  project_dir, "outputs", "hcc_final_model_shap_external",
  "final_model_shap_external_objects.rds"
))
shap5 <- readRDS(file.path(
  project_dir, "outputs", "hcc_final_model_shap_external_16s_only",
  "final_model_shap_external_objects.rds"
))
comparison <- readRDS(file.path(out_dir, "reviewer_16s_only_comparison_objects.rds"))

emb6 <- as.data.frame(shap6$internal_embedding)
emb5 <- as.data.frame(shap5$internal_embedding)
metrics <- as.data.frame(comparison$cluster_metrics)
model_labels <- c("Six cohorts (5 16S + MGX)", "Five 16S cohorts only")

base_theme <- theme_hcc_embedding(base_size = 10) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", color = "#111111", hjust = 0),
    legend.position = "bottom"
  )

make_cluster_plot <- function(dat, title, colors) {
  ggplot(dat, aes(tSNE1, tSNE2, color = Cluster)) +
    geom_point(alpha = 0.70, size = 1.65) +
    scale_color_manual(values = colors, limits = names(colors), drop = FALSE) +
    labs(title = title, x = "t-SNE 1", y = "t-SNE 2", color = "SHAP cluster") +
    base_theme
}

p_clusters <- (
  make_cluster_plot(emb6, model_labels[1], HCC_CLUSTER_COLORS_6) |
    make_cluster_plot(emb5, model_labels[2], HCC_CLUSTER_COLORS_16S)
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "SHAP clusters: four clusters with MGX versus three clusters in 16S-only data",
    subtitle = sprintf(
      "Shared-sample ARI = %.3f; best-overlap agreement = %.1f%%",
      metrics$AdjustedRandIndex,
      100 * metrics$BestOverlapAgreement
    ),
    theme = theme(plot.title = element_text(face = "bold", color = "#111111"))
  ) &
  theme(legend.position = "bottom")


