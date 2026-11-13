#=======================================perdataset Shannon analysis



library(readxl)

file_path <- "revise metadata final liver cancer add metagenomics.xlsx"  

# read Sheet1
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


# perdataset shannon diversity for boxplot
library(vegan)
library(dplyr)
library(tidyr)


otu <- t(rel_abund)

# count Shannon
shannon <- diversity(otu, index = "shannon")

# merge metadata
alpha_df <- data.frame(
  RunID = names(shannon),
  Shannon = shannon
)

alpha_df <- alpha_df %>%
  left_join(metadata, by = "RunID")


stat_res <- alpha_df %>%
  group_by(ProjectID) %>%
  summarise(
    p_value = wilcox.test(Shannon ~ Status)$p.value
  )

library(ggplot2)
library(ggpubr)
p <- ggplot(alpha_df, aes(x = Status, y = Shannon, fill = Status)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "gray30", lwd = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, shape = 21, color = "black") +
  
  scale_fill_manual(values = c("HD" = "#67a8cd", "HCC" = "#f36569")) +
  
  facet_wrap(~ ProjectID, scales = "free") +
  labs(title = "Shannon Diversity between HCC and HD", x = NULL, y = "Shannon diversity") +
  theme_bw() +
  labs(
    y = "Shannon diversity",
    x = ""
  ) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()   
  ) + 
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )
  

p




#=======================================SMD meta-analysis


library(dplyr)
library(tidyr)

meta_input <- alpha_df %>%
  
  group_by(ProjectID, Status) %>%
  summarise(
    n = n(),
    mean = mean(Shannon, na.rm = TRUE),
    sd = sd(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from = Status,
    values_from = c(n, mean, sd)
  ) %>%
  
  rename(
    study = ProjectID,
    
    n_HCC = n_HCC,
    mean_HCC = mean_HCC,
    sd_HCC = sd_HCC,
    
    n_HD = n_HD,
    mean_HD = mean_HD,
    sd_HD = sd_HD
  ) %>%
  
  filter(!is.na(n_HCC) & !is.na(n_HD))


# meta-analysis（SMD）
library(meta)

m <- metacont(
  n.e = n_HCC,
  mean.e = mean_HCC,
  sd.e = sd_HCC,
  
  n.c = n_HD,
  mean.c = mean_HD,
  sd.c = sd_HD,
  
  studlab = study,
  data = meta_input,
  
  sm = "SMD",
  method.tau = "REML",
  
  method.random.ci = "HK"  
)




forest(
  m,
  prediction = TRUE,

  leftcols = c("studlab", "n.e", "mean.e", "sd.e",
               "n.c", "mean.c", "sd.c"),
  
  leftlabs = c("Study",
               "Total", "Mean", "SD",
               "Total", "Mean", "SD"),
  

  rightlabs = c("SMD", "95% CI"),
  

  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.pval.Q = TRUE,
  

  col.diamond = "black",
  col.diamond.lines = "black",
  

  col.square = "grey40",
  col.square.lines = "black",
  

  col.inside = "black",
  

  xlab = "Standardised mean difference",
  

  fontsize = 10,
  

  comb.random = TRUE,
  comb.fixed = F
)






























































#=======================================perdataset Simpson analysis



library(readxl)

file_path <- "revise metadata final liver cancer add metagenomics.xlsx"  # 如果文件在当前工作目录# 或者使用完整路径

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



library(vegan)
library(dplyr)
library(tidyr)


otu <- t(rel_abund)


simpson <- diversity(otu, index = "simpson")   # 修改


alpha_df <- data.frame(
  RunID = names(simpson),
  Simpson = simpson        
)

alpha_df <- alpha_df %>%
  left_join(metadata, by = "RunID")


stat_res <- alpha_df %>%
  group_by(ProjectID) %>%
  summarise(
    p_value = wilcox.test(Simpson ~ Status)$p.value   # 修改
  )

library(ggplot2)
library(ggpubr)
p <- ggplot(alpha_df, aes(x = Status, y = Simpson, fill = Status)) +   
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "gray30", lwd = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, shape = 21, color = "black") +
  
  scale_fill_manual(values = c("HD" = "#67a8cd", "HCC" = "#f36569")) +
  
  facet_wrap(~ ProjectID, scales = "free") +
  labs(title = "Simpson Diversity between HCC and HD", x = NULL, y = "Simpson diversity") +   
  theme_bw() +
  labs(
    y = "Simpson diversity",   
    x = ""
  ) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()   
  ) + 
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )

p



library(dplyr)
library(tidyr)

meta_input <- alpha_df %>%
  

  group_by(ProjectID, Status) %>%
  summarise(
    n = n(),
    mean = mean(Simpson, na.rm = TRUE), 
    sd = sd(Simpson, na.rm = TRUE),      
    .groups = "drop"
  ) %>%
  

  pivot_wider(
    names_from = Status,
    values_from = c(n, mean, sd)
  ) %>%
  

  rename(
    study = ProjectID,
    
    n_HCC = n_HCC,
    mean_HCC = mean_HCC,
    sd_HCC = sd_HCC,
    
    n_HD = n_HD,
    mean_HD = mean_HD,
    sd_HD = sd_HD
  ) %>%
  

  filter(!is.na(n_HCC) & !is.na(n_HD))

library(meta)

m <- metacont(
  n.e = n_HCC,
  mean.e = mean_HCC,
  sd.e = sd_HCC,
  
  n.c = n_HD,
  mean.c = mean_HD,
  sd.c = sd_HD,
  
  studlab = study,
  data = meta_input,
  
  sm = "SMD",
  method.tau = "REML",
  
  method.random.ci = "HK"  
)


forest(
  m,
  prediction = TRUE,
  leftcols = c("studlab", "n.e", "mean.e", "sd.e",
               "n.c", "mean.c", "sd.c"),
  leftlabs = c("Study",
               "Total", "Mean", "SD",
               "Total", "Mean", "SD"),
  rightlabs = c("SMD", "95% CI"),
  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.pval.Q = TRUE,
  col.diamond = "black",
  col.diamond.lines = "black",
  col.square = "grey40",
  col.square.lines = "black",
  col.inside = "black",
  xlab = "Standardised mean difference",
  fontsize = 10,
  comb.random = TRUE,
  comb.fixed = FALSE
)














###：Observed species


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


# perdataset Observed species for boxplot
library(vegan)
library(dplyr)
library(tidyr)


otu <- t(rel_abund)  

# Observed species
observed <- specnumber(otu)   

alpha_df <- data.frame(
  RunID = names(observed),   
  Observed = observed        
)

alpha_df <- alpha_df %>%
  left_join(metadata, by = "RunID")



stat_res <- alpha_df %>%
  group_by(ProjectID) %>%
  summarise(
    p_value = wilcox.test(Observed ~ Status)$p.value   
  )

library(ggplot2)
library(ggpubr)
p <- ggplot(alpha_df, aes(x = Status, y = Observed, fill = Status)) +   
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "gray30", lwd = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, shape = 21, color = "black") +
  
  scale_fill_manual(values = c("HD" = "#67a8cd", "HCC" = "#f36569")) +
  
  facet_wrap(~ ProjectID, scales = "free") +
  labs(title = "Observed species between HCC and HD", x = NULL, y = "Observed species") +   
  theme_bw() +
  labs(
    y = "Observed species",   
    x = ""
  ) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()   
  ) + 
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  )




library(dplyr)
library(tidyr)

meta_input <- alpha_df %>%
  
  
  group_by(ProjectID, Status) %>%
  summarise(
    n = n(),
    mean = mean(Observed, na.rm = TRUE),   
    sd = sd(Observed, na.rm = TRUE),       
    .groups = "drop"
  ) %>%
  

  pivot_wider(
    names_from = Status,
    values_from = c(n, mean, sd)
  ) %>%
  
 
  rename(
    study = ProjectID,
    
    n_HCC = n_HCC,
    mean_HCC = mean_HCC,
    sd_HCC = sd_HCC,
    
    n_HD = n_HD,
    mean_HD = mean_HD,
    sd_HD = sd_HD
  ) %>%
  
  
  filter(!is.na(n_HCC) & !is.na(n_HD))


library(meta)

m <- metacont(
  n.e = n_HCC,
  mean.e = mean_HCC,
  sd.e = sd_HCC,
  
  n.c = n_HD,
  mean.c = mean_HD,
  sd.c = sd_HD,
  
  studlab = study,
  data = meta_input,
  
  sm = "SMD",
  method.tau = "REML",
  
  method.random.ci = "HK"   
)


forest(
  m,
  prediction = TRUE,
  leftcols = c("studlab", "n.e", "mean.e", "sd.e",
               "n.c", "mean.c", "sd.c"),
  leftlabs = c("Study",
               "Total", "Mean", "SD",
               "Total", "Mean", "SD"),
  rightlabs = c("SMD", "95% CI"),
  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.pval.Q = TRUE,
  col.diamond = "black",
  col.diamond.lines = "black",
  col.square = "grey40",
  col.square.lines = "black",
  col.inside = "black",
  xlab = "Standardised mean difference",
  fontsize = 10,
  comb.random = TRUE,
  comb.fixed = FALSE
)

















































