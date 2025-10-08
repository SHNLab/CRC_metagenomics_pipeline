# Title: Final Comprehensive Analysis Pipeline (Blue-White-Red Heatmap)
# Author: abhilasha, aman
# Institute: NABI, Mohali, India
## Date: 2025-09-15
# Description: This is the definitive, all-in-one script for the analysis of
#              flux balance data. It includes a blue-white-red heatmap,
#              PERMANOVA, PCA, and other publication-style figures (1200 dpi).

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these packages installed:
# install.packages(c("tidyverse", "fs", "ggrepel", "pheatmap", "RColorBrewer", "vegan"))

library(tidyverse)
library(fs)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(vegan)      # For PERMANOVA (adonis2)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_FINAL_16spet")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:132)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)
print("✅ Metadata generated.")

# --- Load and Consolidate Metabolite Flux Data ---
all_scfa_files <- fs::dir_ls(path = diet_path, glob = "**/SCFA_Harmful_Metabolite_Profile*.csv", recurse = TRUE)
if (length(all_scfa_files) == 0) stop("CRITICAL ERROR: No SCFA/Metabolite files found.")

all_metabolite_data <- all_scfa_files %>%
  purrr::map_dfr(~ {
    file_path <- .x
    file_name <- fs::path_file(file_path)
    diet_name <- fs::path_split(file_path)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE) %>%
      select(starts_with("Metabolite_"), starts_with("WeightedFlux_")) %>%
      pivot_longer(everything(),
                   names_to = c(".value", "id"),
                   names_pattern = "(.+)_(\\d+)") %>%
      filter(!is.na(Metabolite) & !is.na(WeightedFlux)) %>%
      mutate(
        sample_id = sample_id,
        diet = diet_name,
        Metabolite = str_to_title(str_replace_all(Metabolite, "_", " "))
      ) %>%
      select(sample_id, diet, Metabolite, WeightedFlux)
  })

# --- Final Merged Data with Unique IDs ---
merged_flux_data <- all_metabolite_data %>%
  inner_join(metadata, by = "sample_id") %>%
  mutate(sample_diet_id = paste(sample_id, diet, sep = "_")) # Unique ID for each row
print("✅ Metabolite flux data loaded and merged.")


# 4. STATISTICAL ANALYSIS
# ----------------------------------------------------
# --- Differential Metabolite Analysis ---
stats_results <- merged_flux_data %>%
  group_by(diet, Metabolite) %>%
  filter(n_distinct(group) == 2) %>%
  summarise(
    p_value = wilcox.test(WeightedFlux ~ group)$p.value,
    mean_healthy = mean(WeightedFlux[group == "Healthy"], na.rm = TRUE),
    mean_crc = mean(WeightedFlux[group == "CRC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_healthy > 0 & mean_crc > 0) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    log2_fold_change = log2(mean_crc / mean_healthy)
  )
readr::write_csv(stats_results, fs::path(output_dir, "metabolite_statistics_summary.csv"))
print("✅ Differential metabolite analysis complete.")

# --- Multivariate Analysis Prep ---
multivariate_data_wide <- merged_flux_data %>%
  pivot_wider(id_cols = c(sample_diet_id, diet, group), names_from = Metabolite, values_from = WeightedFlux, values_fill = 0)

metabolite_matrix <- multivariate_data_wide %>% select(-sample_diet_id, -diet, -group)
sample_info <- multivariate_data_wide %>% select(sample_diet_id, diet, group)

# Remove zero-variance columns (critical fix for PCA and PERMANOVA)
variances <- apply(metabolite_matrix, 2, var, na.rm = TRUE)
metabolite_matrix_filtered <- metabolite_matrix[, variances > 0]

# --- PERMANOVA Statistical Test ---
print("--- Running PERMANOVA with Euclidean Distance ---")
permanova_result <- adonis2(metabolite_matrix_filtered ~ group * diet, data = sample_info, permutations = 999, method = "euclidean")
permanova_table <- as.data.frame(permanova_result)
print(permanova_table)
write.csv(permanova_table, fs::path(output_dir, "permanova_results_detailed.csv"))


# 5. FIGURE AND TABLE GENERATION
# ----------------------------------------------------

# --- Figure 1: Heatmap of Metabolite Fluxes (Blue-White-Red colors) ---
annotation_col <- data.frame(Condition = sample_info$group, Diet = sample_info$diet)
rownames(annotation_col) <- sample_info$sample_diet_id

heatmap_matrix <- t(scale(metabolite_matrix_filtered))
colnames(heatmap_matrix) <- sample_info$sample_diet_id

png(fs::path(output_dir, "Figure1_Metabolite_Heatmap.png"), width = 12, height = 10, units = "in", res = 1200)
pheatmap(
  mat = heatmap_matrix,
  main = "Global Metabolite Flux Profiles (Z-score)",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  # ** CHANGED COLOR PALETTE to Blue-White-Red **
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  fontsize = 14,
  cutree_rows = 2
)
dev.off()
print("✅ Figure 1 (Heatmap) saved.")


# --- Figure 2: Volcano Plots for Each Diet ---
volcano_dir <- fs::path(output_dir, "Figure2_Volcano_Plots")
dir_create(volcano_dir)

for(current_diet in unique(stats_results$diet)) {
  plot_data <- stats_results %>%
    filter(diet == current_diet) %>%
    mutate(
      significance = case_when(p_adj < 0.05 & abs(log2_fold_change) > 1 ~ "Significant", TRUE ~ "Not Significant"),
      label = if_else(significance == "Significant", Metabolite, "")
    )
  volcano_plot <- ggplot(plot_data, aes(x = log2_fold_change, y = -log10(p_adj))) +
    geom_point(aes(color = significance), alpha = 0.8, size = 3) +
    geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 15) +
    scale_color_manual(values = c("Significant" = "#F7374F", "Not Significant" = "grey")) +
    labs(title = paste("Dysregulated Metabolites in", current_diet, "Diet"), subtitle = "CRC vs. Healthy", x = "Log2 Fold Change (CRC / Healthy)", y = "-log10 (Adjusted P-value)") +
    theme_bw(base_size = 14) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted")
  ggsave(filename = fs::path(volcano_dir, paste0("Volcano_", current_diet, ".png")), plot = volcano_plot, width = 8, height = 7, dpi = 1200)
}
print("✅ Figure 2 (Volcano Plots) saved.")


# --- Figure 3: Diet Intervention "Rescue" Plot ---
rescue_metabolites <- c("Butyrate", "Hydrogen Sulfide")
rescue_data <- merged_flux_data %>%
  filter(Metabolite %in% rescue_metabolites) %>%
  group_by(diet, group, Metabolite) %>%
  summarise(mean_flux = mean(WeightedFlux, na.rm = TRUE), se = sd(WeightedFlux, na.rm = TRUE) / sqrt(n()), .groups = "drop")

rescue_plot <- ggplot(rescue_data, aes(x = diet, y = mean_flux, color = group, group = group)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean_flux - se, ymax = mean_flux + se), width = 0.2) +
  facet_wrap(~Metabolite, scales = "free_y", labeller = as_labeller(c(Butyrate = "Beneficial: Butyrate Production", `Hydrogen Sulfide` = "Harmful: Hydrogen Sulfide Production"))) +
  scale_color_manual(values = c("CRC" = "#F7374F", "Healthy" = "#229799")) +
  labs(title = "Evaluating Diet's Ability to 'Rescue' CRC Metabolic Profile", subtitle = "Goal: Diets that bring the red line (CRC) closer to the blue line (Healthy)", x = "Diet Intervention", y = "Mean Weighted Flux") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
ggsave(filename = fs::path(output_dir, "Figure3_Diet_Rescue_Effect.png"), plot = rescue_plot, width = 12, height = 7, dpi = 1200)
print("✅ Figure 3 (Rescue Plot) saved.")


# --- Figure 4: PCA of Metabolite Profiles ---
pca_fit <- prcomp(metabolite_matrix_filtered, scale. = TRUE, center = TRUE)
pca_plot_data <- as_tibble(pca_fit$x) %>%
  bind_cols(sample_info)

pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group, shape = diet)) +
  geom_point(size = 5, alpha = 0.8) +
  labs(
    title = "PCA of Metabolite Profiles Across Diets",
    subtitle = "Samples cluster based on overall metabolic similarity",
    x = paste0("PC1 (", round(summary(pca_fit)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_fit)$importance[2,2]*100, 1), "%)"),
    color = "Health Status", shape = "Diet"
  ) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("CRC" = "#F7374F", "Healthy" = "#229799")) +
  coord_fixed()
ggsave(filename = fs::path(output_dir, "Figure4_PCA_Metabolite_Profiles.png"), plot = pca_plot, width = 12, height = 9, dpi = 1200)
print("✅ Figure 4 (PCA Plot) saved.")


# --- Final Summary Table ---
rescue_summary_table <- rescue_data %>%
  select(-se) %>%
  pivot_wider(names_from = group, values_from = mean_flux) %>%
  mutate(flux_difference = Healthy - CRC)
readr::write_csv(rescue_summary_table, fs::path(output_dir, "rescue_effect_summary_table.csv"))
print("✅ Rescue effect summary table saved.")

# --- Script End ---
print("--- 🎉 Analysis complete! All files saved in 'publication_style_analysis_FINAL' folder. ---")