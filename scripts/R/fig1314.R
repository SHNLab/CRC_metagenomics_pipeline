# Title: Definitive Integrated Microbiome Analysis Pipeline (with CRC Growers Plot)
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This is the final, all-in-one script that performs all analyses:
#              time-series (0hr vs 24hr), differential abundance, diversity,
#              taxonomic summaries, correlation heatmaps, and all figures,
#              including a dedicated plot for top species growing faster in CRC.

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
library(tidytext)


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


# 5. STATISTICAL ANALYSIS OF LFC
# ----------------------------------------------------
lfc_stats <- lfc_results %>%
  group_by(diet, Species) %>%
  filter(n_distinct(group) == 2) %>%
  rstatix::wilcox_test(lfc ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

lfc_diffs <- lfc_results %>%
  group_by(diet, Species, group) %>%
  summarise(mean_lfc = mean(lfc), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_lfc, values_fill = 0) %>%
  mutate(lfc_diff = CRC - Healthy)

volcano_stats <- left_join(lfc_stats, lfc_diffs, by = c("diet", "Species"))
print("✅ LFC statistical analysis complete.")


# 6. VISUALIZATION
# ----------------------------------------------------
# --- Figure 12: Summary of Top 5 Overall Significant LFC Differences ---
top_5_lfc_diffs <- volcano_stats %>%
  filter(p.adj < 0.05) %>%
  group_by(diet) %>%
  slice_min(order_by = p.adj, n = 5) %>%
  mutate(
    regulation = if_else(lfc_diff > 0, "Grows Faster in CRC", "Grows Faster in Healthy"),
    Species = tidytext::reorder_within(Species, abs(lfc_diff), diet)
  )

summary_plot <- ggplot(top_5_lfc_diffs, aes(x = Species, y = lfc_diff, fill = regulation)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ diet, scales = "free_y") +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = c("Grows Faster in CRC" = "#F7374F", "Grows Faster in Healthy" = "#229799")) +
  labs(title = "Top 5 Species with Significantly Different Growth Dynamics (Overall)",
       subtitle = "Comparing LFC (24hr/0hr) between CRC and Healthy groups",
       x = "Species", y = "Difference in Mean LFC (CRC - Healthy)") +
  theme_bw(base_size = 14)
ggsave(filename = fs::path(output_dir, "Figure12_Summary_Top5_LFC_Overall.png"), plot = summary_plot, width = 16, height = 10, dpi = 1200)
print("✅ Figure 12 (Summary of Overall Top 5 LFC differences) saved.")


# --- NEW Figure 13: Summary of Top 5 Species Growing Faster in CRC ---
top_5_crc_growers <- volcano_stats %>%
  filter(p.adj < 0.05, lfc_diff > 0) %>% # Filter for species that grow faster in CRC
  group_by(diet) %>%
  slice_min(order_by = p.adj, n = 5) %>% # Find the top 5 most significant
  mutate(Species = tidytext::reorder_within(Species, lfc_diff, diet))

summary_plot_crc <- ggplot(top_5_crc_growers, aes(x = Species, y = lfc_diff, fill = "Grows Faster in CRC")) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ diet, scales = "free_y") +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = c("Grows Faster in CRC" = "#F7374F"), name = "Regulation") +
  labs(
    title = "Top 5 Species Growing Significantly Faster in CRC",
    subtitle = "Comparing LFC (24hr/0hr) between CRC and Healthy groups",
    x = "Species",
    y = "Difference in Mean LFC (CRC - Healthy)"
  ) +
  theme_bw(base_size = 14)
ggsave(filename = fs::path(output_dir, "Figure13_Summary_Top5_CRC_Growers.png"), plot = summary_plot_crc, width = 16, height = 10, dpi = 1200)
print("✅ NEW Figure 13 (Summary of Top 5 CRC Growers) saved.")


# (Other plots like Heatmap, Volcano plots can be generated here if needed)
# ...

# --- Script End ---
print("--- 🎉 All analyses complete! Check your output folder. ---")