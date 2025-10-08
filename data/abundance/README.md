# Abundance Data

Description:
------------
This folder contains species- and genus-level abundance tables generated from Kraken2 and Bracken outputs. These tables provide normalized microbial composition for individual samples and combined cohort-level summaries.  

Folder Structure:
-----------------
abundance/
  CRC/
    merged_genus_output.txt     # Combined genus-level abundance across all CRC patients
    merged_species_output.txt   # Combined species-level abundance across all CRC patients
    merged_kraken_report.txt    # Kraken2 summary report for all CRC patients

  Healthy/
    merged_genus_output.txt     # Combined genus-level abundance across all healthy samples
    merged_species_output.txt   # Combined species-level abundance across all healthy samples
    merged_kraken_report.txt    # Kraken2 summary report for all healthy samples
   

File Description:
-----------------
- `merged_genus_output.txt`: Combined genus-level abundance table across all samples in the group.
- `merged_species_output.txt`: Combined species-level abundance table across all samples in the group.
- `merged_kraken_report.txt`: Kraken2 summary report combining all samples in the group.
- `New Text Document.txt`: Optional text file containing metadata, notes, or readme information.

Usage:
------
These abundance tables are used for:
1. Microbial community analyses (alpha/beta diversity)
2. Personalized metabolic modeling (FBA simulations)
