# Title: Definitive Personalized and Collective Diet Prediction Analysis
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-17
# Description: This definitive script includes a robustness check for missing
#              samples, generates personalized plots for ALL CRC samples with the
#              Set2 color palette, and creates a collective summary heatmap
#              with statistical validation.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these libraries:
# install.packages(c("tidyverse", "fs", "readxl", "patchwork", "RColorBrewer", "rstatix", "pheatmap"))

library(tidyverse)
library(fs)
library(readxl)
library(patchwork)    # For combining plots into a single figure
library(RColorBrewer)
library(rstatix)
library(pheatmap)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
baseline_path <- fs::path(base_path, "brakcen_filtered")
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_PREDICTIONS")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
# Using the corrected, continuous list from your script
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)
print("✅ Metadata generated.")

# --- Load 0hr Baseline Abundance Data ---
baseline_files <- fs::dir_ls(baseline_path, regexp = "\\.xlsx$")
if (length(baseline_files) == 0) stop("CRITICAL ERROR: No .xlsx files found.")
baseline_data <- baseline_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    sample_id <- stringr::str_remove(file_name, "\\.xlsx$")
    read_excel(.x) %>%
      rename(Model = 1, RelAbundance = 2) %>%
      mutate(
        sample_id = sample_id,
        Species = stringr::str_replace(Model, "\\.mat$", "") %>% stringr::str_replace_all("_", " ")
      )
  })
print("✅ Baseline (0hr) abundance data loaded.")

# --- Load 24hr Simulated Metabolite Data ---
all_scfa_files <- fs::dir_ls(path = diet_path, glob = "**/SCFA_Harmful_Metabolite_Profile*.csv", recurse = TRUE)
if (length(all_scfa_files) == 0) stop("CRITICAL ERROR: No SCFA/Metabolite files found.")
all_metabolite_data <- all_scfa_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    diet_name <- fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      select(starts_with("Metabolite_"), starts_with("WeightedFlux_")) %>%
      pivot_longer(everything(), names_to = c(".value", "id"), names_pattern = "(.+)_(\\d+)") %>%
      filter(!is.na(Metabolite) & !is.na(WeightedFlux)) %>%
      mutate(
        sample_id = sample_id,
        diet = diet_name,
        Metabolite = str_to_title(str_replace_all(Metabolite, "_", " "))
      ) %>%
      select(sample_id, diet, Metabolite, WeightedFlux)
  })
merged_flux_data <- all_metabolite_data %>%
  inner_join(metadata, by = "sample_id")
print("✅ Simulated (24hr) metabolite data loaded.")


# 4. GENERATE INDIVIDUAL PREDICTION PLOTS FOR ALL CRC SAMPLES
# ----------------------------------------------------
# --- Establish the "Healthy Target" Profile ---
healthy_target_profile <- merged_flux_data %>%
  filter(group == "Healthy") %>%
  group_by(diet, Metabolite) %>%
  summarise(healthy_avg_flux = mean(WeightedFlux, na.rm = TRUE), .groups = "drop")

# Create a subdirectory for the individual plots
individual_plots_dir <- fs::path(output_dir, "Individual_Patient_Reports")
if (!fs::dir_exists(individual_plots_dir)) dir_create(individual_plots_dir)

# --- Define color palette ---
set2_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:length(unique(merged_flux_data$diet))]
names(set2_colors) <- unique(merged_flux_data$diet)

# --- Loop through each CRC sample ---
for (current_sample in crc_ids) {
  
  cat(paste("\n--- Processing Sample:", current_sample, "---\n"))
  
  # **FIX**: Check if the sample exists in the loaded data before proceeding
  if (!(current_sample %in% merged_flux_data$sample_id)) {
    warning(paste("Sample", current_sample, "not found in the data. Skipping report generation."))
    next # Skip to the next iteration of the loop
  }
  
  sample_flux_data <- merged_flux_data %>%
    filter(sample_id == current_sample) %>%
    left_join(healthy_target_profile, by = c("diet", "Metabolite")) %>%
    mutate(
      response = case_when(
        Metabolite %in% c("Butyrate", "Propionate", "Acetate") & WeightedFlux > healthy_avg_flux ~ "Rescued",
        Metabolite %in% c("Hydrogen Sulfide", "Indole", "Ammonia") & WeightedFlux < healthy_avg_flux ~ "Rescued",
        TRUE ~ "Not Rescued"
      )
    )
  
  plot_a <- ggplot(sample_flux_data, aes(x = diet, y = WeightedFlux, fill = diet)) +
    geom_col(alpha = 1) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    geom_hline(aes(yintercept = healthy_avg_flux), linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = set2_colors) +
    labs(subtitle = "A) Predicted Metabolite Flux (Dashed line = Healthy Avg.)", x = NULL, y = "Flux") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  plot_b <- ggplot(sample_flux_data, aes(x = diet, y = Metabolite, fill = response)) +
    geom_tile(color = "white", lwd = 1.5) +
    scale_fill_manual(values = c("Rescued" = "purple", "Not Rescued" = "grey90")) +
    labs(subtitle = "B) Rescue Status", x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_c <- baseline_data %>%
    filter(sample_id == current_sample) %>%
    slice_max(order_by = RelAbundance, n = 10) %>%
    ggplot(aes(x = reorder(Species, RelAbundance), y = RelAbundance)) +
    geom_col(fill = "steelblue", alpha = 1) +
    coord_flip() +
    labs(subtitle = "C) Top 10 Abundant Species (Baseline)", x = NULL, y = "Relative Abundance") +
    theme_bw(base_size = 14)
  
  final_plot <- (plot_a | (plot_b / plot_c)) +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(title = paste("Personalized Dietary Prediction for Sample:", current_sample))
  
  ggsave(
    filename = fs::path(individual_plots_dir, paste0("Prediction_Report_", current_sample, ".png")),
    plot = final_plot, width = 18, height = 10, dpi = 1200
  )
  print(paste("✅ Report saved for", current_sample))
}


# 5. GENERATE COLLECTIVE SUMMARY PLOT
# ----------------------------------------------------
# --- Calculate Average "Rescue Effect" ---
collective_data <- merged_flux_data %>%
  left_join(healthy_target_profile, by = c("diet", "Metabolite")) %>%
  mutate(lfc_vs_healthy = log2((WeightedFlux + 1e-9) / (healthy_avg_flux + 1e-9)))

# --- Perform Statistical Test ---
collective_stats <- collective_data %>%
  filter(group == "CRC") %>%
  group_by(diet, Metabolite) %>%
  rstatix::wilcox_test(lfc_vs_healthy ~ 1, mu = 0) %>%
  rstatix::adjust_pvalue(method = "fdr")

# --- Prepare data for Heatmap ---
heatmap_data <- collective_data %>%
  filter(group == "CRC") %>%
  group_by(diet, Metabolite) %>%
  summarise(mean_lfc = mean(lfc_vs_healthy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = diet, values_from = mean_lfc) %>%
  column_to_rownames("Metabolite")

# Prepare significance matrix
significance_matrix <- collective_stats %>%
  mutate(stars = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))) %>%
  select(diet, Metabolite, stars) %>%
  pivot_wider(names_from = diet, values_from = stars, values_fill = "") %>%
  column_to_rownames("Metabolite")
significance_matrix <- significance_matrix[rownames(heatmap_data), colnames(heatmap_data)]


# --- Generate Heatmap ---
png(fs::path(output_dir, "Figure14_Collective_Rescue_Effect_Heatmap.png"), width = 12, height = 10, units = "in", res = 1200)
pheatmap(
  mat = heatmap_data,
  main = "Collective Diet Rescue Effect on CRC Metabolite Profiles",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  display_numbers = significance_matrix,
  fontsize = 14,
  fontsize_number = 20,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("Lower in CRC\n(Rescued)", "-1", "Healthy Avg.", "+1", "Higher in CRC\n(Not Rescued)")
)
dev.off()
print("✅ NEW Figure 14 (Collective Rescue Heatmap) saved.")


# --- Script End ---
print("\n--- 🎉 All personalized and collective analyses are complete! ---")