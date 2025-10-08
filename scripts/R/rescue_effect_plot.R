# Title: Diet Rescue Effect Plot for All Metabolites
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-15
# Description: This script loads the metabolite flux data and generates a
#              separate "rescue" plot for every metabolite, visualizing the
#              difference between CRC and Healthy groups across all diets.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(ggpubr)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
# Create a new, dedicated folder for these plots
output_dir <- fs::path(base_path, "Metabolite_Rescue_Plots")
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

# --- Load Metabolite Data ---
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
print("✅ Metabolite flux data loaded.")


# 4. CALCULATE SUMMARY DATA AND GENERATE PLOTS
# ----------------------------------------------------
# --- Calculate mean and standard error for all metabolites ---
rescue_data_all <- merged_flux_data %>%
  group_by(diet, group, Metabolite) %>%
  summarise(
    mean_flux = mean(WeightedFlux, na.rm = TRUE),
    se = sd(WeightedFlux, na.rm = TRUE) / sqrt(n()), # Standard Error
    .groups = "drop"
  )

# --- Get the list of all unique metabolites ---
all_metabolites <- unique(rescue_data_all$Metabolite)

# --- Loop through each metabolite and create a separate plot ---
for (current_metabolite in all_metabolites) {
  
  cat(paste("Generating plot for:", current_metabolite, "\n"))
  
  # Filter the data for the current metabolite
  plot_data <- rescue_data_all %>%
    filter(Metabolite == current_metabolite)
  
  # Create the plot
  rescue_plot <- ggplot(plot_data, aes(x = diet, y = mean_flux, color = group, group = group)) +
    geom_point(size = 4) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = mean_flux - se, ymax = mean_flux + se), width = 0.2) +
    scale_color_manual(values = c("CRC" = "#F7374F", "Healthy" = "#229799")) +
    labs(
      title = paste("Diet's Effect on", current_metabolite, "Production"),
      subtitle = "Comparing Mean Flux in CRC vs. Healthy Groups",
      x = "Diet Intervention",
      y = "Mean Weighted Flux"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot with a unique name
  ggsave(
    filename = fs::path(output_dir, paste0("Rescue_Plot_", current_metabolite, ".png")),
    plot = rescue_plot,
    width = 10,
    height = 7,
    dpi = 1200
  )
}

# --- Script End ---
print(paste("\n--- 🎉 All metabolite plots are saved in:", output_dir, "---"))