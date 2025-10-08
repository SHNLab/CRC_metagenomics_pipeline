# Title: Final CRC Rescue Status Plot & Summary Generator
# Author: abhilasha
# Institute: NABI, Mohali, India
# Description: This script generates a single figure showing the rescue status
#              for every individual CRC patient, with the x-axis correctly sorted
#              numerically. It also exports the raw data to a CSV and a
#              summarized analysis to a TXT file.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)
library(readxl)

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

# --- Define the specific list of metabolites to include ---
target_metabolites <- c(
  "Butyrate", "Propionate", "Acetate", "Lactate", "Formate", "Succinate",
  "Valerate", "Indole", "Phenol", "Putrescine", "Ammonia", "Trimethylamine", "Hydrogen Sulfide"
)

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
print("\n--- Analyzing rescue status and preparing outputs... ---")

# --- Establish the "Healthy Target" Profile ---
healthy_target_profile <- merged_flux_data %>%
  filter(group == "Healthy") %>%
  group_by(diet, Metabolite) %>%
  summarise(healthy_avg_flux = mean(WeightedFlux, na.rm = TRUE), .groups = "drop")

# --- Calculate Rescue Status and correctly order patients ---
crc_rescue_data <- merged_flux_data %>%
  filter(group == "CRC", Metabolite %in% target_metabolites) %>%
  filter(!is.na(Metabolite)) %>%
  left_join(healthy_target_profile, by = c("diet", "Metabolite")) %>%
  mutate(
    response = case_when(
      Metabolite %in% c("Butyrate", "Propionate", "Acetate", "Lactate", "Formate", "Succinate", "Valerate") & WeightedFlux > healthy_avg_flux ~ "Rescued",
      Metabolite %in% c("Hydrogen Sulfide", "Indole", "Ammonia", "Trimethylamine", "Phenol", "Putrescine") & WeightedFlux < healthy_avg_flux ~ "Rescued",
      TRUE ~ "Not Rescued"
    ),
    Metabolite = factor(Metabolite, levels = target_metabolites),
    # **NEW**: Create a numeric column for sorting
    patient_num = as.numeric(str_extract(sample_id, "[0-9]+"))
  ) %>%
  # **NEW**: Arrange by the patient number
  arrange(patient_num) %>%
  # **NEW**: Set the factor levels for sample_id to the sorted order
  mutate(sample_id = factor(sample_id, levels = unique(.data$sample_id)))


# --- Save the full underlying data as a CSV file ---
write_csv(crc_rescue_data, fs::path(output_dir, "Combined_CRC_Rescue_Status_Data.csv"))
print(paste("✅ Full rescue status data saved to CSV."))


# --- **NEW**: Generate and save the summary TXT file ---
# Summary 1: Count of rescued patients per metabolite for each diet
summary_by_metabolite <- crc_rescue_data %>%
  filter(response == "Rescued") %>%
  group_by(diet, Metabolite) %>%
  summarise(patients_rescued_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = diet, values_from = patients_rescued_count, values_fill = 0)

# Summary 2: "Rescue Score" - total rescued metabolites for each patient per diet
summary_by_patient <- crc_rescue_data %>%
  filter(response == "Rescued") %>%
  group_by(diet, sample_id) %>%
  summarise(rescued_metabolites_score = n(), .groups = "drop") %>%
  pivot_wider(names_from = diet, values_from = rescued_metabolites_score, values_fill = 0)

# Write summaries to a text file
summary_file_path <- fs::path(output_dir, "Rescue_Summary_Results.txt")
sink(summary_file_path) # Redirects output to the file

cat("====================================================\n")
cat("         RESCUE STATUS SUMMARY RESULTS\n")
cat("====================================================\n\n")

cat("--- Summary 1: Patients Rescued per Metabolite (out of 28) ---\n")
cat("This table shows how many patients showed a 'Rescued' status for each metabolite under each diet.\n\n")
print(summary_by_metabolite, n = Inf) # n=Inf ensures all rows are printed

cat("\n\n--- Summary 2: Rescue Score per Patient ---\n")
cat("This table shows the total number of metabolites that were 'Rescued' for each patient by each diet.\n\n")
print(summary_by_patient, n = Inf)

sink() # Closes the connection, stopping the redirect

print(paste("✅ Summary results saved to:", summary_file_path))


# --- Create the combined plot with sorted X-axis ---
combined_rescue_plot <- ggplot(crc_rescue_data, aes(x = sample_id, y = Metabolite, fill = response)) +
  geom_tile(color = "white", lwd = 0.5) +
  facet_wrap(~ diet, ncol = 3, scales = "free_x") +
  scale_fill_manual(
    name = "Status",
    values = c("Rescued" = "#687FE5", "Not Rescued" = "#DCDCDC"),
    labels = c("Rescued", "Not Rescued")
  ) +
  labs(
    title = "Collective Rescue Status of Metabolites Across All CRC Patients",
    subtitle = "Each column represents an individual patient, faceted by dietary intervention.",
    x = "CRC Patient ID (Sorted)",
    y = "Metabolite"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

# --- Save the plot at 1200 DPI ---
ggsave(
  filename = fs::path(output_dir, "Figure_Combined_CRC_Rescue_Status_Sorted_1200dpi.png"),
  plot = combined_rescue_plot,
  width = 18,
  height = 12,
  dpi = 1200
)

print(paste("✅ Sorted rescue plot saved successfully!"))