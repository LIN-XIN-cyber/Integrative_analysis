

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
local_lib <- Sys.getenv("HCC_R_LIB", unset = file.path(project_dir, "R_library"))
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
source(file.path(project_dir, "analysis", "plot_style.R"))

source_rds <- if (length(args) >= 1) args[1] else file.path(project_dir, "outputs", "hcc_shap_rank", "cohort_shap_rank_objects.rds")
out_dir <- if (length(args) >= 2) args[2] else dirname(source_rds)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(source_rds)
rank_table <- as.data.frame(obj$mean_abs_shap_ranks)
instability <- as.data.frame(obj$rank_instability)
cohorts <- sort(unique(rank_table$Cohort))
pipeline_names <- c("Before MMUPHin", "After MMUPHin")

plot_table <- rank_table %>%
  left_join(instability %>% select(Pipeline, Feature, PanelSize, InstabilityRank),
    by = c("Pipeline", "Feature")) %>%
  mutate(
    DisplayRank = if_else(Present, as.numeric(Rank), as.numeric(PanelSize + 1L)),
    Cohort = factor(Cohort, levels = cohorts)
  )

genus_palette <- hcc_genus_rank_colors(sort(unique(plot_table$Feature)))

make_rank_plot <- function(dat, pipeline) {
  dat <- filter(dat, Pipeline == pipeline)
  features <- sort(unique(dat$Feature))
  panel_size <- unique(dat$PanelSize)
  pal <- genus_palette[features]
  top_drivers <- dat %>% distinct(Feature, InstabilityRank) %>%
    filter(InstabilityRank <= 5) %>% pull(Feature)

  ggplot(dat, aes(Cohort, DisplayRank, group = Feature, color = Feature)) +
    geom_line(color = "white", linewidth = 2.0, alpha = 0.92) +
    geom_line(aes(linewidth = InstabilityRank <= 5), alpha = 0.94) +
    geom_point(
      data = filter(dat, Present), aes(fill = Feature),
      shape = 21, color = "#1F2937", stroke = 0.45, size = 3.2, alpha = 1,
      show.legend = FALSE
    ) +
    geom_point(
      data = filter(dat, !Present), shape = 4, size = 4.0, stroke = 1.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, labels = sub("^g__", "", names(pal))) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_linewidth_manual(values = c(`FALSE` = 0.85, `TRUE` = 1.55), guide = "none") +
    scale_y_reverse(
      breaks = seq_len(panel_size + 1L),
      labels = c(as.character(seq_len(panel_size)), "Absent"),
      minor_breaks = NULL
    ) +
    labs(
      title = paste0(pipeline, ": genus importance-rank stability across cohorts"),
      subtitle = NULL,
      x = "Cohort", y = "Mean |SHAP| rank (1 = most important)",
      color = "Genus",
      caption = paste0(
        "Within-cohort RF mean |SHAP| rank; thicker lines mark the five most unstable taxa.\n",
        "Five most unstable genera: ", paste(sub("^g__", "", top_drivers), collapse = ", "), "\n",
        "× = genus absent from that cohort (zero prevalence); absent taxa were excluded from that cohort's RF and SHAP calculation."
      )
    ) +
    theme_hcc_publication(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(
        face = "bold", color = "#111111", size = 15, hjust = 0,
        margin = margin(b = 5)
      ),
      plot.caption = element_text(
        color = "#555555", hjust = 0, size = 8.5, lineheight = 1.13,
        margin = margin(t = 6)
      ),
      legend.position = "right", legend.text = element_text(size = 8.3),
      legend.key.height = grid::unit(0.52, "cm"), legend.key.width = grid::unit(0.9, "cm")
    ) +
    guides(color = guide_legend(ncol = ifelse(length(features) > 10, 2, 1), byrow = TRUE,
      override.aes = list(linewidth = 2.0, alpha = 1)))
}

save_rank_plot <- function(plot, slug) {
  png_path <- file.path(fig_dir, paste0("cohort_shap_rank_", slug, ".png"))
  pdf_path <- file.path(fig_dir, paste0("cohort_shap_rank_", slug, ".pdf"))
  editable_pdf_path <- hcc_editable_pdf_path(pdf_path)
  temp_png <- file.path(fig_dir, paste0("cohort_shap_rank_", slug, "_palette_tmp.png"))
  temp_pdf <- file.path(fig_dir, paste0("cohort_shap_rank_", slug, "_palette_tmp.pdf"))
  temp_editable_pdf <- hcc_editable_pdf_path(temp_pdf)

  ggsave(temp_png, plot, width = 16.5, height = 9.5, dpi = 320, bg = "white")
  ggsave(temp_pdf, plot, width = 16.5, height = 9.5, device = grDevices::cairo_pdf, bg = "white")
  ggsave(temp_pdf, plot, width = 16.5, height = 9.5, device = hcc_editable_pdf_device, bg = "white")

  replacements <- c(
    file.copy(temp_png, png_path, overwrite = TRUE),
    file.copy(temp_pdf, pdf_path, overwrite = TRUE),
    file.copy(temp_editable_pdf, editable_pdf_path, overwrite = TRUE)
  )
  if (!all(replacements)) stop("Could not replace one or more SHAP-rank figure files.")
  unlink(c(temp_png, temp_pdf, temp_editable_pdf))
}

for (pipeline in pipeline_names) {
  slug <- ifelse(pipeline == "Before MMUPHin", "before_mmuphin", "after_mmuphin")
  p <- make_rank_plot(plot_table, pipeline)
  save_rank_plot(p, slug)
}

color_key <- data.frame(
  Feature = names(genus_palette),
  Genus = sub("^g__", "", names(genus_palette)),
  Color = unname(genus_palette),
  stringsAsFactors = FALSE
)
utils::write.csv(
  color_key,
  file.path(out_dir, "tables_csv", "genus_color_key.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Replotted existing SHAP-rank objects without refitting models: ", out_dir)
