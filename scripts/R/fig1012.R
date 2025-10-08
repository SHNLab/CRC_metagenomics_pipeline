# Title: Definitive Integrated Microbiome Analysis Pipeline (Final Version)
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This is the final, all-in-one script that performs all analyses:
#              time-series (0hr vs 24hr), differential abundance, diversity,
#              taxonomic summaries, correlation heatmaps, and all figures.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these packages installed:
# install.packages(c("tidyverse", "fs", "readxl", "ggpubr", "vegan", "pheatmap", "RColorBrewer", "rstatix", "ggrepel", "tidytext"))

library(tidyverse)
library(fs)
library(readxl)
library(ggpubr)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(rstatix)
library(ggrepel)
library(tidytext) # <-- ADDED THIS LIBRARY FOR reorder_within()


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
baseline_path <- fs::path(base_path, "brakcen_filtered")
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

# --- Load 0hr Baseline Data ---
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
        time = "0hr",
        Species = stringr::str_replace(Model, "\\.mat$", "") %>% stringr::str_replace_all("_", " ")
      )
  })
print("✅ Baseline (0hr) data loaded.")

# --- Load 24hr Simulated Data ---
simulated_files <- fs::dir_ls(path = diet_path, glob = "**/updated_abundances*.csv", recurse = TRUE)
if (length(simulated_files) == 0) stop("CRITICAL ERROR: No 24hr abundance files found.")
simulated_data <- simulated_files %>%
  purrr::map_dfr(~ {
    file_name <- fs::path_file(.x)
    diet_name <- fs::path_split(.x)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(file_name, "(di_p|hi_p)[0-9]+")
    if(is.na(sample_id)) return(NULL)
    readr::read_csv(.x, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        diet = diet_name,
        sample_id = sample_id,
        time = "24hr",
        Species = stringr::str_replace(Model, "\\.mat$", "") %>% stringr::str_replace_all("_", " ")
      )
  })
print("✅ Simulated (24hr) data loaded.")

# --- Combine Data ---
baseline_replicated <- tidyr::crossing(baseline_data, diet = unique(simulated_data$diet))
time_series_data <- bind_rows(baseline_replicated, simulated_data)
print("✅ Combined 0hr and 24hr data.")


# 4. LOG FOLD CHANGE (LFC) CALCULATION
# ----------------------------------------------------
lfc_data_wide <- time_series_data %>%
  pivot_wider(id_cols = c(sample_id, diet, Species), names_from = time, values_from = RelAbundance, values_fill = 0) %>%
  filter(`0hr` > 0, `24hr` > 0)

lfc_results <- lfc_data_wide %>%
  mutate(lfc = log2((`24hr` + 1e-6) / (`0hr` + 1e-6))) %>%
  inner_join(metadata, by = "sample_id")

readr::write_csv(lfc_results, fs::path(output_dir, "species_lfc_results.csv"))
print("✅ Log Fold Change calculated.")


# 5. VISUALIZATION AND ANALYSIS OF LFC
# ----------------------------------------------------
# --- Figure 10: Heatmap of Log Fold Change ---
top_dynamic_species <- lfc_results %>%
  group_by(Species) %>%
  summarise(mean_abs_lfc = mean(abs(lfc))) %>%
  slice_max(order_by = mean_abs_lfc, n = 50) %>%
  pull(Species)

heatmap_matrix_lfc <- lfc_results %>%
  filter(Species %in% top_dynamic_species) %>%
  mutate(sample_diet_id = paste(sample_id, diet, sep = "_")) %>%
  pivot_wider(id_cols = Species, names_from = sample_diet_id, values_from = lfc, values_fill = 0) %>%
  column_to_rownames("Species")

annotation_col <- lfc_results %>%
  select(sample_id, diet, group) %>% distinct() %>%
  mutate(sample_diet_id = paste(sample_id, diet, sep = "_")) %>%
  select(sample_diet_id, diet, group) %>%
  column_to_rownames("sample_diet_id")

png(fs::path(output_dir, "Figure10_LFC_Heatmap.png"), width = 16, height = 12, units = "in", res = 1200)
pheatmap(mat = heatmap_matrix_lfc, main = "Log Fold Change (24hr vs 0hr) of Top 50 Dynamic Species",
         annotation_col = annotation_col, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         fontsize = 12, fontsize_row = 8, cluster_cols = TRUE)
dev.off()
print("✅ Figure 10 (LFC Heatmap) saved.")


# --- Figure 11: Volcano Plots for ALL Diets (Loop) ---
lfc_stats <- lfc_results %>%
  group_by(diet, Species) %>%
  filter(n_distinct(group) == 2) %>%
  rstatix::wilcox_test(lfc ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

lfc_diffs <- lfc_results %>%
  group_by(diet, Species, group) %>%
  summarise(mean_lfc = mean(lfc), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_lfc) %>%
  mutate(lfc_diff = CRC - Healthy)

volcano_stats <- left_join(lfc_stats, lfc_diffs, by = c("diet", "Species"))

volcano_dir <- fs::path(output_dir, "Figure11_Volcano_Plots_All_Diets")
dir_create(volcano_dir)

for(current_diet in unique(volcano_stats$diet)) {
  plot_data <- volcano_stats %>%
    filter(diet == current_diet) %>%
    mutate(
      significance = if_else(p.adj < 0.05, "Significant", "Not Significant"),
      label = if_else(significance == "Significant", Species, "")
    )
  volcano_lfc_plot <- ggplot(plot_data, aes(x = lfc_diff, y = -log10(p.adj))) +
    geom_point(aes(color = significance), size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = label), max.overlaps = 15) +
    scale_color_manual(values = c("Significant" = "#FF894F", "Not Significant" = "#FDEBD0")) +
    labs(title = paste("Species with Different Growth Dynamics in", current_diet, "Diet"),
         subtitle = "Comparing LFC (24hr/0hr) between CRC and Healthy",
         x = "Difference in Mean LFC (CRC - Healthy)") +
    theme_bw(base_size = 14)
  ggsave(filename = fs::path(volcano_dir, paste0("Volcano_LFC_", current_diet, ".png")),
         plot = volcano_lfc_plot, width = 10, height = 8, dpi = 1200)
}
print("✅ Figure 11 (Volcano plots for all diets) saved.")


# --- Figure 12: Summary of Top 5 Significant LFC Differences ---
top_5_lfc_diffs <- volcano_stats %>%
  filter(p.adj < 0.05) %>%
  group_by(diet) %>%
  slice_min(order_by = p.adj, n = 5) %>%
  mutate(
    regulation = if_else(lfc_diff > 0, "Grows Faster in CRC", "Grows Faster in Healthy"),
    Species = tidytext::reorder_within(Species, -p.adj, diet) # Correctly call the function
  )

summary_plot <- ggplot(top_5_lfc_diffs, aes(x = Species, y = lfc_diff, fill = regulation)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ diet, scales = "free_y") +
  tidytext::scale_x_reordered() + # Use the partner function for plotting
  scale_fill_manual(values = c("Grows Faster in CRC" = "#F7374F", "Grows Faster in Healthy" = "#229799")) +
  labs(
    title = "Top 5 Species with Significantly Different Growth Dynamics",
    subtitle = "Comparing LFC (24hr/0hr) between CRC and Healthy groups",
    x = "Species",
    y = "Difference in Mean LFC (CRC - Healthy)"
  ) +
  theme_bw(base_size = 14)

ggsave(filename = fs::path(output_dir, "Figure12_Summary_Top5_LFC_Differences.png"),
       plot = summary_plot, width = 16, height = 10, dpi = 1200)
print("✅ Figure 12 (Summary plot of Top 5 LFC differences) saved.")


# --- Script End ---
print("--- 🎉 All analyses complete! Check your output folder. ---")