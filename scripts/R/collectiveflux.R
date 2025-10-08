# Title: Final CRC Flux Plot & Summary Generator
# Institute: NABI, Mohali, India
# Description: This script generates a single figure showing metabolite flux
#              for every individual CRC patient, with the x-axis correctly sorted
#              numerically. It also exports the raw data to a CSV and a
#              summarized analysis to a TXT file.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(readxl)
library(RColorBrewer)

# 2. CONFIGURATION: Define Paths & Create Output Directory
# ----------------------------------------------------
# !!! IMPORTANT: Please update this path to your working directory !!!
base_path <- "/Volumes/Backups of /bacarena"

# --- Define sub-directories ---
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_PREDICTIONS")

# --- Create the output directory if it doesn't exist ---
if (!fs::dir_exists(output_dir)) {
  dir_create(output_dir)
  print(paste("Output directory created at:", output_dir))
}

# 3. DATA PREPARATION
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)
print("✅ Metadata generated.")

# --- Load and Process Simulated Metabolite Data ---
all_scfa_files <- fs::dir_ls(path = diet_path, glob = "**/SCFA_Harmful_Metabolite_Profile*.csv", recurse = TRUE)

if (length(all_scfa_files) == 0) stop("CRITICAL ERROR: No metabolite files found in the specified 'diet_path'.")

all_metabolite_data <- all_scfa_files %>%
  purrr::map_dfr(~ {
    file_path <- .x
    diet_name <- fs::path_split(file_path)[[1]] %>% purrr::pluck(length(.) - 2)
    sample_id <- stringr::str_extract(fs::path_file(file_path), "(di_p|hi_p)[0-9]+")
    
    if(is.na(sample_id)) return(NULL)
    
    readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE) %>%
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

print("✅ Simulated metabolite data loaded and processed.")


# 4. ANALYSIS, PLOTTING, AND EXPORT
# ----------------------------------------------------
print("\n--- Analyzing flux data and preparing outputs... ---")

# --- Establish the "Healthy Target" Profile ---
healthy_target_profile <- merged_flux_data %>%
  filter(group == "Healthy") %>%
  group_by(diet, Metabolite) %>%
  summarise(healthy_avg_flux = mean(WeightedFlux, na.rm = TRUE), .groups = "drop")

# --- Define color palette ---
set2_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[1:length(unique(merged_flux_data$diet))]
names(set2_colors) <- unique(merged_flux_data$diet)

# --- Filter data for CRC samples and correctly order patients ---
crc_flux_data <- merged_flux_data %>%
  filter(group == "CRC") %>%
  left_join(healthy_target_profile, by = c("diet", "Metabolite")) %>%
  # **NEW**: Create a numeric column for sorting
  mutate(patient_num = as.numeric(str_extract(sample_id, "[0-9]+"))) %>%
  # **NEW**: Arrange by the patient number
  arrange(patient_num) %>%
  # **NEW**: Set the factor levels for sample_id to the sorted order
  mutate(sample_id = factor(sample_id, levels = unique(.data$sample_id)))

# --- Save the full underlying data as a CSV file ---
write_csv(crc_flux_data, fs::path(output_dir, "Combined_CRC_Flux_Data.csv"))
print(paste("✅ Full flux data saved to CSV."))


# --- **NEW**: Generate and save the summary TXT file ---
# Summary 1: Average flux per metabolite for each diet
summary_by_metabolite <- crc_flux_data %>%
  group_by(diet, Metabolite) %>%
  summarise(mean_flux = mean(WeightedFlux, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = diet, values_from = mean_flux)

# Summary 2: Average flux per patient for each diet
summary_by_patient <- crc_flux_data %>%
  group_by(diet, sample_id) %>%
  summarise(mean_flux_across_metabolites = mean(WeightedFlux, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = diet, values_from = mean_flux_across_metabolites)

# Write summaries to a text file
summary_file_path <- fs::path(output_dir, "Flux_Summary_Results.txt")
sink(summary_file_path) # Redirects output to the file

cat("====================================================\n")
cat("            METABOLITE FLUX SUMMARY RESULTS\n")
cat("====================================================\n\n")

cat("--- Summary 1: Mean Flux per Metabolite ---\n")
cat("This table shows the average flux for each metabolite across all 28 patients, for each diet.\n\n")
print(summary_by_metabolite, n = Inf) # n=Inf ensures all rows are printed

cat("\n\n--- Summary 2: Mean Flux per Patient ---\n")
cat("This table shows the average flux across all metabolites for each patient, for each diet.\n\n")
print(summary_by_patient, n = Inf)

sink() # Closes the connection, stopping the redirect

print(paste("✅ Summary results saved to:", summary_file_path))


# --- Create the combined plot with sorted X-axis ---
combined_crc_plot <- ggplot(crc_flux_data, aes(x = sample_id, y = WeightedFlux)) +
  geom_col(aes(fill = diet), show.legend = FALSE) +
  geom_hline(aes(yintercept = healthy_avg_flux), linetype = "dashed", color = "red", linewidth = 1) +
  facet_grid(Metabolite ~ diet, scales = "free_y") +
  scale_fill_manual(values = set2_colors) +
  labs(
    title = "Individual Metabolite Flux Across All CRC Patients by Diet",
    subtitle = "Each bar is one patient. Red dashed line indicates the Healthy cohort average.",
    x = "CRC Patient ID (Sorted)",
    y = "Predicted Metabolite Flux"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

# --- Save the plot at 1200 DPI ---
ggsave(
  filename = fs::path(output_dir, "Figure_Combined_CRC_Metabolite_Flux_Sorted_1200dpi.png"),
  plot = combined_crc_plot,
  width = 20,
  height = 16,
  dpi = 1200
)




print(paste("✅ Sorted flux plot saved successfully!"))

#END#