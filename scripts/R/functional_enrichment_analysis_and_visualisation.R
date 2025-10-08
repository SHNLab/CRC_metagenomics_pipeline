# Title: Functional Category Enrichment Analysis and Visualization
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-19
# Description: This self-contained script loads raw data, maps metabolites to
#              health-related functional categories, and visualizes which functions
#              are enriched in CRC vs. Healthy using a publication-quality bubble plot.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
# Make sure you have these packages installed:
# install.packages(c("tidyverse", "fs", "rstatix", "viridis"))

library(tidyverse)
library(fs)
library(rstatix)
library(viridis)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
base_path <- "/Volumes/Backups of /bacarena"
diet_path <- fs::path(base_path, "diet")
output_dir <- fs::path(base_path, "publication_style_analysis_FUNCTIONAL_BUBBLEPLOT")
if (!fs::dir_exists(output_dir)) dir_create(output_dir)


# 3. DEFINE FUNCTIONAL METABOLITE CATEGORIES
# ----------------------------------------------------
# This is where we map each metabolite to its known biological function.
# This list can be customized based on literature and your specific interests.
functional_categories <- list(
  `SCFA Production` = c("Butyrate", "Propionate", "Acetate", "Valerate"),
  `Inflammation-Associated` = c("Hydrogen Sulfide", "Putrescine", "Ammonia", "Indole"),
  `Energy Metabolism (TCA Cycle)` = c("Succinate", "Formate"),
  `Amino Acid Metabolism` = c("Trimethylamine", "Indole", "Putrescine", "Ammonia"),
  `Lactic Acid Fermentation` = c("Lactate")
)
print("✅ Functional metabolite categories defined.")


# 4. LOAD AND PREPARE DATA (Self-Contained)
# ----------------------------------------------------
# --- Generate Metadata ---
crc_ids <- paste0("di_p", c(1:28))
healthy_ids <- paste0("hi_p", 29:133)
metadata <- bind_rows(
  tibble(sample_id = crc_ids, group = "CRC"),
  tibble(sample_id = healthy_ids, group = "Healthy")
)

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


# 5. PERFORM STATISTICAL ANALYSIS
# ----------------------------------------------------
# --- Perform Wilcoxon test to find significant metabolites ---
metabolite_stats <- merged_flux_data %>%
  group_by(diet, Metabolite) %>%
  filter(n_distinct(group) == 2 &&
           (var(WeightedFlux[group == "CRC"], na.rm = TRUE) > 0 || var(WeightedFlux[group == "Healthy"], na.rm = TRUE) > 0)) %>%
  rstatix::wilcox_test(WeightedFlux ~ group) %>%
  rstatix::adjust_pvalue(method = "fdr")

# --- Calculate Log2 Fold Change ---
mean_fluxes <- merged_flux_data %>%
  group_by(diet, Metabolite) %>%
  summarise(
    mean_healthy = mean(WeightedFlux[group == "Healthy"], na.rm = TRUE),
    mean_crc = mean(WeightedFlux[group == "CRC"], na.rm = TRUE),
    .groups = "drop"
  )

# --- Combine stats with fold change and map to functional categories ---
analysis_data <- metabolite_stats %>%
  left_join(mean_fluxes, by = c("diet", "Metabolite")) %>%
  mutate(
    log2fc = log2((mean_crc + 1e-9) / (mean_healthy + 1e-9)),
    Category = setNames(rep(names(functional_categories), lengths(functional_categories)), unlist(functional_categories))[Metabolite]
  ) %>%
  filter(!is.na(Category))
print("✅ Data mapped to functional categories.")


# 6. CALCULATE ENRICHMENT AND VISUALIZE
# ----------------------------------------------------
# --- Calculate Category-Level Statistics for the plot ---
category_summary <- analysis_data %>%
  group_by(diet, Category) %>%
  summarise(
    num_significant = sum(p.adj < 0.05, na.rm = TRUE),
    avg_log2fc = mean(log2fc, na.rm = TRUE),
    total_metabolites = n()
  ) %>%
  ungroup() %>%
  mutate(proportion_significant = num_significant / total_metabolites)

write_csv(category_summary, fs::path(output_dir, "functional_category_summary.csv"))
print("✅ Functional category summary calculated.")


# --- Generate the Bubble Plot ---
functional_bubble_plot <- ggplot(category_summary,
                                 aes(x = diet, y = Category, size = proportion_significant, color = avg_log2fc)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "#3B9AB2", mid = "white", high = "#F21A00", midpoint = 0,
                        name = "Avg. Log2 Fold Change\n(CRC vs. Healthy)") +
  scale_size(range = c(4, 20), name = "Proportion of\nSignificant Metabolites") +
  labs(
    title = "Functional Enrichment of Metabolic Categories by Diet",
    subtitle = "Red indicates enrichment in CRC; Blue indicates enrichment in Healthy",
    x = "Diet",
    y = "Functional Category"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted")
  )

ggsave(filename = fs::path(output_dir, "Figure_Functional_Category_Enrichment.png"),
       plot = functional_bubble_plot, width = 14, height = 8, dpi = 1200)

print("✅ Functional enrichment bubble plot saved.")

# --- Script End ---