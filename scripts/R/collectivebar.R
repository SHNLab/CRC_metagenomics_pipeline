# Title: Collective CRC Group Metabolite Profile Analysis
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-17
# Description: This script generates summary bar plots for each metabolite,
#              showing the mean and standard deviation of flux for the entire
#              CRC cohort across all six diets, compared to the healthy average.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(RColorBrewer)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_COLLECTIVE_SUMMARY")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. METADATA AND DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)
print("✅ Metadata generated.")

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


# 4. CALCULATE GROUP AVERAGES AND GENERATE PLOTS
# ----------------------------------------------------
# --- Create a dedicated folder for the plots ---
collective_plots_dir <- fs::path(output_dir, "Collective_Metabolite_Plots")
if (!fs::dir_exists(collective_plots_dir)) dir_create(collective_plots_dir)

# --- Calculate the average flux for the Healthy group (our target baseline) ---
healthy_target_profile <- merged_flux_data %>%
  filter(group == "Healthy") %>%
  group_by(diet, Metabolite) %>%
  summarise(healthy_avg_flux = mean(WeightedFlux, na.rm = TRUE), .groups = "drop")

# --- Calculate the mean and SD for the CRC group ---
crc_summary <- merged_flux_data %>%
  filter(group == "CRC") %>%
  group_by(diet, Metabolite) %>%
  summarise(
    mean_flux_crc = mean(WeightedFlux, na.rm = TRUE),
    sd_flux_crc = sd(WeightedFlux, na.rm = TRUE),
    .groups = "drop"
  )

# --- Combine CRC summary with the Healthy target for plotting ---
plot_data_full <- crc_summary %>%
  left_join(healthy_target_profile, by = c("diet", "Metabolite"))

# --- Define color palette ---
set2_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:length(unique(plot_data_full$diet))]
names(set2_colors) <- unique(plot_data_full$diet)

# --- Initialize the summary report ---
report_lines <- list()
report_lines[[1]] <- "--- Collective Metabolite Summary Report ---"

# --- Loop through each metabolite to create a separate plot and report ---
for (current_metabolite in unique(plot_data_full$Metabolite)) {
  
  cat(paste("Generating plot and report for:", current_metabolite, "\n"))
  
  plot_data_metabolite <- plot_data_full %>%
    filter(Metabolite == current_metabolite)
  
  # --- Add to Summary Report ---
  report_lines[[length(report_lines) + 1]] <- paste0("\n## Metabolite: ", current_metabolite)
  report_lines[[length(report_lines) + 1]] <- "  **Comparison of Mean Flux (CRC vs. Healthy Average):**"
  for(i in 1:nrow(plot_data_metabolite)) {
    row <- plot_data_metabolite[i,]
    report_lines[[length(report_lines) + 1]] <- paste0(
      "    - ", row$diet, " Diet: CRC Mean = ", round(row$mean_flux_crc, 2),
      " (SD = ", round(row$sd_flux_crc, 2), "), Healthy Mean = ", round(row$healthy_avg_flux, 2)
    )
  }
  
  # --- Generate the Bar Plot ---
  metabolite_plot <- ggplot(plot_data_metabolite, aes(x = diet, y = mean_flux_crc, fill = diet)) +
    geom_col(alpha = 1) +
    geom_errorbar(aes(ymin = mean_flux_crc - sd_flux_crc, ymax = mean_flux_crc + sd_flux_crc),
                  width = 0.25, linewidth = 0.8) +
    geom_hline(aes(yintercept = healthy_avg_flux), linetype = "dashed", color = "red", linewidth = 1.2) +
    facet_wrap(~ diet, scales = "free_x") + # Facet to give each bar its own space
    scale_fill_manual(values = set2_colors) +
    labs(
      title = paste("CRC Cohort: Mean Production of", current_metabolite),
      subtitle = "Bars represent Mean ± SD for the CRC group across all diets. Red dashed line is the Healthy cohort average.",
      x = "Diet",
      y = "Mean Weighted Flux"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(), # Remove redundant x-axis labels
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 14)
    )
  
  ggsave(
    filename = fs::path(collective_plots_dir, paste0("Collective_Plot_", current_metabolite, ".png")),
    plot = metabolite_plot,
    width = 12,
    height = 8,
    dpi = 1200
  )
}


# 5. WRITE THE FINAL SUMMARY REPORT
# ----------------------------------------------------
output_file_path <- fs::path(output_dir, "Collective_Metabolite_Summary.txt")
writeLines(unlist(report_lines), con = output_file_path)

print(paste("\n--- 🎉 Success! Your summary report has been saved to:", output_file_path, "---"))
print(paste("--- All individual metabolite plots are saved in:", collective_plots_dir, "---"))

# --- Script End ---