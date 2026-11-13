

library(readxl)

file_path <- "revise metadata final liver cancer add metagenomics.xlsx"

metadata <- read_excel(
  path = file_path,
  sheet = 1
)

library(readr)

tax_profile <- read_tsv("all_six_genus_level_tables.tsv")
tax_profile <- as.data.frame(tax_profile)
names(tax_profile)[1] <- "taxon"
rownames(tax_profile) <- tax_profile$taxon

tax_profile <- tax_profile[,-1]

tax_profile <- tax_profile[,names(tax_profile) %in% metadata$RunID]

tax_profile <- tax_profile[rowSums(tax_profile != 0) > 0, ]

rel_abund <- tax_profile

prev <- rowMeans(rel_abund > 0)
rel_abund_filt <- rel_abund[prev > 0.1, ]

mean_abund <- rowMeans(rel_abund_filt)
rel_abund_filt <- rel_abund_filt[mean_abund > 1e-4, ]

rel_abund_filt <- sweep(rel_abund_filt, 2, colSums(rel_abund_filt), "/")

abundance <- as.data.frame(t(rel_abund_filt))

abundance <- as.data.frame(t(abundance))

metadata$ProjectID <- as.factor(metadata$ProjectID)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$RunID

samples_meta <- rownames(metadata)
abundance <- abundance[, samples_meta]

library(readr)
library(dplyr)
library(Maaslin2)
library(ggplot2)
library(vegan)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(readxl)



batch_cor_filtered_abundance <- read_tsv("batch_cor_genus_profile.tsv")

batch_cor_filtered_abundance <- as.data.frame(batch_cor_filtered_abundance)

rownames(batch_cor_filtered_abundance) <- batch_cor_filtered_abundance$...1
batch_cor_filtered_abundance <- batch_cor_filtered_abundance[,-1]
batch_cor_filtered_abundance <- as.data.frame(t(batch_cor_filtered_abundance))

batch_cor_filtered_abundance <- as.data.frame(t(batch_cor_filtered_abundance))

original_abundance <- as.data.frame(abundance)

metadata$Status <- as.factor(metadata$Status)
metadata$RunID <- as.character(metadata$RunID)
metadata$ProjectID <- as.factor(metadata$ProjectID)


original_abundance_maaslin = Maaslin2(
  input_data = original_abundance, 
  input_metadata = metadata, 
  normalization = "NONE",
  transform = "LOG",
  output = "original_abundance" , 
  fixed_effects = c( "Status"), 
  reference = c("Status,HCC"),
  plot_heatmap = TRUE
)

batch_cor_filtered_abundance_data = Maaslin2(
  input_data = batch_cor_filtered_abundance, 
  input_metadata = metadata, 
  normalization = "NONE",
  transform = "LOG",
  output = "batch_cor_filtered_abundance" , 
  fixed_effects = c( "Status"), 
  reference = c("Status,HCC"),
  plot_heatmap = TRUE
)



maaslin_all_transform <- read_tsv("original_abundance/transformed_maaslin.tsv")

maaslin_all_batch <- read_tsv("batch_cor_filtered_abundance/batch_maaslin.tsv")

maaslin_sig_transform <- read_tsv("original_abundance/transformed_maaslin_significant_results.tsv")

maaslin_sig_batch <- read_tsv("batch_cor_filtered_abundance/batch_maaslin_significant_results.tsv")

maaslin_all_transform$feature <- gsub("g__","", maaslin_all_transform$feature)
maaslin_all_batch$feature <- gsub("g__","", maaslin_all_batch$feature)
maaslin_sig_transform$feature <- gsub("g__","", maaslin_sig_transform$feature)
maaslin_sig_batch$feature <- gsub("g__","", maaslin_sig_batch$feature)

library(ggplot2)
library(dplyr)
library(ggrepel)

maaslin_all_transform <- maaslin_all_transform %>%
  mutate(
    log10qval = -log10(qval),  
    direction = ifelse(coef > 0, "Enriched in HD", "Enriched in HCC"),
    abs_coef = abs(coef),
    color_group = ifelse(qval > 0.05, "nonsig",
                         ifelse(coef > 0, "Enriched in HD", "Enriched in HCC"))
  )

top_features <- maaslin_all_transform %>%
  filter(qval <= 0.05) %>%
  arrange(desc(abs_coef)) %>%
  slice_head(n = 10)

volcano_plot_transform <- ggplot(maaslin_all_transform, aes(x = coef, y = log10qval)) +
  geom_point(
    data = subset(maaslin_all_transform, qval > 0.05),
    aes(size = abs_coef),
    color = "grey75", alpha = 0.4
  ) +
  geom_point(
    data = subset(maaslin_all_transform, qval <= 0.05),
    aes(color = direction, size = abs_coef),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c("Enriched in HCC" = "#FAC1C3",
               "Enriched in HD"  = "#C2DCEB")
  ) +
  scale_size_continuous(range = c(2, 6)) +
  geom_text_repel(
    data = top_features,
    aes(label = feature),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Volcano Plot of Maaslin2 Results",
    subtitle = "original relab log transformed",
    x = "Coefficient",
    y = "-log10(q-value)",
    color = "Direction",
    size = "|Coefficient|"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title = element_text(face = "bold")
  ) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

maaslin_all_batch <- maaslin_all_batch %>%
  mutate(
    log10qval = -log10(qval),  
    direction = ifelse(coef > 0, "Enriched in HD", "Enriched in HCC"),
    abs_coef = abs(coef),
    color_group = ifelse(qval > 0.05, "nonsig",
                         ifelse(coef > 0, "Enriched in HD", "Enriched in HCC"))
  )

top_features <- maaslin_all_batch %>%
  filter(qval <= 0.05) %>%
  arrange(desc(abs_coef)) %>%
  slice_head(n = 10)

volcano_plot_batch <- ggplot(maaslin_all_batch, aes(x = coef, y = log10qval)) +
  geom_point(
    data = subset(maaslin_all_batch, qval > 0.05),
    aes(size = abs_coef),
    color = "grey75", alpha = 0.4
  ) +
  geom_point(
    data = subset(maaslin_all_batch, qval <= 0.05),
    aes(color = direction, size = abs_coef),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c("Enriched in HCC" = "#FAC1C3",
               "Enriched in HD"  = "#C2DCEB")
  ) +
  scale_size_continuous(range = c(2, 6)) +
  geom_text_repel(
    data = top_features,
    aes(label = feature),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Volcano Plot of Maaslin2 Results",
    subtitle = "batch correction",
    x = "Coefficient",
    y = "-log10(q-value)",
    color = "Direction",
    size = "|Coefficient|"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title = element_text(face = "bold")
  ) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

