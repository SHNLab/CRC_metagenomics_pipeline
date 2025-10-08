# Title: Definitive Results Summary Generator
# Author: abhilasha
# Institute: NABI, Mohali, India
# Date: 2025-09-17
# Description: This definitive script reads all previously generated analysis
#              files from the specified location, extracts key statistical
#              values, and compiles them into a single, formatted text report.

# 1. SETUP: Load necessary libraries
# ----------------------------------------------------
library(tidyverse)
library(fs)


# 2. CONFIGURATION: Define Paths
# ----------------------------------------------------
# This should be the main project directory
base_path <- "/Volumes/Backups of /bacarena"

# **UPDATED**: Path to the main folder containing all your results
main_results_dir <- fs::path(base_path, "publication_style_analysis_FINAL_16spet")


# 3. INITIALIZE THE REPORT
# ----------------------------------------------------
report <- list()
report[[length(report) + 1]] <- "==========================================================="
report[[length(report) + 1]] <- " COMPREHENSIVE ANALYSIS REPORT: CRC vs. HEALTHY MICROBIOME"
report[[length(report) + 1]] <- "==========================================================="
report[[length(report) + 1]] <- paste("Report generated on:", Sys.time())


# 4. OVERALL COMMUNITY ANALYSIS (PERMANOVA)
# ----------------------------------------------------
report[[length(report) + 1]] <- "\n\n## Section 1: Overall Metabolic Profile Differences"
report[[length(report) + 1]] <- "-----------------------------------------------------------"

permanova_file <- fs::path(main_results_dir, "permanova_results_detailed.csv")
if (fs::file_exists(permanova_file)) {
  # **FIX**: Read the CSV and correctly name the first column 'Term'
  permanova_data <- read_csv(permanova_file, show_col_types = FALSE) %>%
    rename(Term = ...1)
  
  group_effect <- permanova_data %>% filter(Term == "group")
  diet_effect <- permanova_data %>% filter(Term == "diet")
  interaction_effect <- permanova_data %>% filter(Term == "group:diet")
  
  report[[length(report) + 1]] <- "\n### PERMANOVA (Metabolic Profile):"
  report[[length(report) + 1]] <- paste0("- Health Status (CRC vs. Healthy) explains ", round(group_effect$R2 * 100, 1), "% of the variance (p = ", group_effect$`Pr(>F)`, ").")
  report[[length(report) + 1]] <- paste0("- Diet explains ", round(diet_effect$R2 * 100, 1), "% of the variance (p = ", diet_effect$`Pr(>F)`, ").")
  report[[length(report) + 1]] <- paste0("- The interaction between diet and health status explains ", round(interaction_effect$R2 * 100, 1), "% of the variance (p = ", interaction_effect$`Pr(>F)`, ").")
  report[[length(report) + 1]] <- "**Interpretation**: Both health status and diet are highly significant drivers of the overall metabolic output."
} else {
  report[[length(report) + 1]] <- "\n- PERMANOVA results file not found."
}


# 5. DIET-BY-DIET ANALYSIS OF ALL METABOLITES
# ----------------------------------------------------
report[[length(report) + 1]] <- "\n\n## Section 2: Diet-Specific Metabolite Profile Analysis"
report[[length(report) + 1]] <- "-----------------------------------------------------------"
metabolite_stats_file <- fs::path(main_results_dir, "metabolite_statistics_summary.csv")
if (fs::file_exists(metabolite_stats_file)) {
  metabolite_stats <- read_csv(metabolite_stats_file, show_col_types = FALSE)
  all_diets <- unique(metabolite_stats$diet)
  
  for (d in all_diets) {
    report[[length(report) + 1]] <- paste0("\n### Diet: ", d)
    diet_data <- metabolite_stats %>% filter(diet == d)
    
    sig_metabolites <- diet_data %>% filter(p_adj < 0.05)
    
    # Calculate log2 fold change using the mean values
    diet_data <- diet_data %>%
      mutate(log2fc = log2((mean_crc + 1e-9) / (mean_healthy + 1e-9)))
    
    if(nrow(sig_metabolites) > 0) {
      report[[length(report) + 1]] <- "  **Significant Metabolite Changes (CRC vs. Healthy):**"
      for(i in 1:nrow(sig_metabolites)) {
        row <- sig_metabolites[i,]
        lfc_row <- diet_data %>% filter(Metabolite == row$Metabolite)
        direction <- if_else(lfc_row$log2fc > 0, "HIGHER in CRC", "LOWER in CRC")
        report[[length(report) + 1]] <- paste0("    - ", row$Metabolite, ": ", direction, " (log2FC = ", round(lfc_row$log2fc, 2), ", p.adj = ", format.pval(row$p_adj, digits = 2), ")")
      }
    } else {
      report[[length(report) + 1]] <- "  - No metabolites were significantly different between CRC and Healthy for this diet."
    }
  }
} else {
  report[[length(report) + 1]] <- "\n- Metabolite statistics file not found."
}


# 6. TIME-SERIES ANALYSIS (0hr vs 24hr)
# ----------------------------------------------------
report[[length(report) + 1]] <- "\n\n## Section 3: Species Growth Dynamics (24hr vs. 0hr)"
report[[length(report) + 1]] <- "-----------------------------------------------------------"
timeseries_stats_file <- fs::path(main_results_dir, "species_lfc_results.csv")
if (fs::file_exists(timeseries_stats_file)) {
  lfc_stats <- read_csv(timeseries_stats_file, show_col_types = FALSE)
  all_diets_ts <- unique(lfc_stats$diet)
  
  for (d in all_diets_ts) {
    report[[length(report) + 1]] <- paste0("\n### Diet: ", d)
    diet_lfc_data <- lfc_stats %>% filter(diet == d)
    
    lfc_diffs <- diet_lfc_data %>%
      group_by(Species, group) %>%
      summarise(mean_lfc = mean(lfc), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean_lfc, values_fill = 0) %>%
      mutate(lfc_diff = CRC - Healthy)
    
    top_crc_growers <- lfc_diffs %>% filter(lfc_diff > 0) %>% arrange(desc(lfc_diff)) %>% head(3)
    if(nrow(top_crc_growers) > 0) {
      report[[length(report) + 1]] <- "  **Top Species Growing Faster in CRC:**"
      for(i in 1:nrow(top_crc_growers)) {
        report[[length(report) + 1]] <- paste0("    - ", top_crc_growers$Species[i], " (LFC difference = +", round(top_crc_growers$lfc_diff[i], 2), ")")
      }
    }
    
    top_healthy_growers <- lfc_diffs %>% filter(lfc_diff < 0) %>% arrange(lfc_diff) %>% head(3)
    if(nrow(top_healthy_growers) > 0) {
      report[[length(report) + 1]] <- "  **Top Species Growing Faster in Healthy:**"
      for(i in 1:nrow(top_healthy_growers)) {
        report[[length(report) + 1]] <- paste0("    - ", top_healthy_growers$Species[i], " (LFC difference = ", round(top_healthy_growers$lfc_diff[i], 2), ")")
      }
    }
  }
} else {
  report[[length(report) + 1]] <- "\n- Time-series statistics file not found."
}


# 7. FUNCTIONAL ENRICHMENT SUMMARY
# ----------------------------------------------------
report[[length(report) + 1]] <- "\n\n## Section 4: Functional Enrichment Summary"
report[[length(report) + 1]] <- "-----------------------------------------------------------"
functional_file <- fs::path(main_results_dir, "functional_enrichment_results.csv")
if(fs::file_exists(functional_file)){
  gsea_results <- read_csv(functional_file, show_col_types = FALSE)
  report[[length(report) + 1]] <- "**Interpretation**: A positive Normalized Enrichment Score (NES) indicates the function is enriched in the CRC group. A negative score indicates enrichment in the Healthy group."
  
  for(d in unique(gsea_results$diet)) {
    report[[length(report) + 1]] <- paste0("\n### Diet: ", d)
    diet_gsea <- gsea_results %>% filter(diet == d)
    for(i in 1:nrow(diet_gsea)){
      row <- diet_gsea[i,]
      significance <- if_else(row$padj < 0.05, paste0("(Significant, p.adj = ", format.pval(row$padj, digits=2), ")"), "(Not Significant)")
      report[[length(report) + 1]] <- paste0("  - ", row$pathway, ": NES = ", round(row$NES, 2), " ", significance)
    }
  }
} else {
  report[[length(report) + 1]] <- "\n- Functional enrichment results file not found."
}


# 8. WRITE THE FINAL REPORT
# ----------------------------------------------------
output_file_path <- fs::path(base_path, "Comprehensive_Results_Summary.txt")
writeLines(unlist(report), con = output_file_path)

print(paste("--- 🎉 Success! Your comprehensive report has been saved to:", output_file_path, "---"))