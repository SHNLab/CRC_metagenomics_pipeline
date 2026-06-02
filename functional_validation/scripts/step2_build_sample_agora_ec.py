#!/usr/bin/env python3
"""
Step 2: Build per-sample AGORA EC sets.

For each usable sample (26 CRC + 103 Healthy = 129):
  - Read patient .xlsx to get list of AGORA models for that sample
  - Look up each model's ECs from the parsed agora_model_ecs.tsv
  - Union all ECs -> "AGORA EC set for this sample"

Output: /scratch/amans/concern3_analysis/per_sample/sample_agora_ecs.tsv
Format: sample_id, ec_number  (long format, one row per sample-EC pair)
"""

import os
import sys
import pandas as pd
from collections import defaultdict

# Paths
AGORA_ECS_TSV = "/scratch/amans/concern3_analysis/parsing/agora_model_ecs.tsv"
MAPPING_TSV = "/scratch/amans/bsim/sample_mapping_tab.tsv"
PATIENT_DIR = "/scratch/amans/bsim/patient"
OUTPUT_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_agora_ecs.tsv"
LOG_FILE = "/scratch/amans/concern3_analysis/logs/step2_build_sample_agora_ec.log"

def log(msg, also_print=True):
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')
    if also_print:
        print(msg)

def main():
    open(LOG_FILE, 'w').close()
    log("=== Step 2: per-sample AGORA EC sets ===")
    
    # Load AGORA model->EC mapping
    log(f"Loading {AGORA_ECS_TSV}")
    df_agora = pd.read_csv(AGORA_ECS_TSV, sep='\t')
    log(f"  {len(df_agora)} model-EC rows, {df_agora['model_name'].nunique()} unique models, {df_agora['ec_number'].nunique()} unique ECs")
    
    # Build dict: model_name -> set of ECs
    model_to_ecs = defaultdict(set)
    for _, row in df_agora.iterrows():
        model_to_ecs[row['model_name']].add(row['ec_number'])
    
    # Load sample mapping
    log(f"Loading {MAPPING_TSV}")
    df_map = pd.read_csv(MAPPING_TSV, sep='\t')
    log(f"  {len(df_map)} mapping entries")
    
    # Build sample -> xlsx_notation lookup
    sample_to_xlsx = dict(zip(df_map['sample_id'], df_map['xlxs_notation']))
    sample_to_group = dict(zip(df_map['sample_id'], df_map['group']))
    
    # Process each sample with available xlsx
    rows_out = []
    samples_processed = 0
    samples_skipped = []
    
    for sample_id, xlsx_notation in sample_to_xlsx.items():
        xlsx_path = os.path.join(PATIENT_DIR, f"{xlsx_notation}.xlsx")
        if not os.path.exists(xlsx_path):
            samples_skipped.append((sample_id, xlsx_notation, "no_xlsx"))
            continue
        
        try:
            df_pat = pd.read_excel(xlsx_path)
        except Exception as e:
            samples_skipped.append((sample_id, xlsx_notation, f"read_error: {str(e)[:80]}"))
            continue
        
        # First column is model filenames (e.g., "Acidaminococcus_fermentans_DSM_20731.xml")
        model_files = df_pat.iloc[:, 0].astype(str).tolist()
        
        # Union ECs across this sample's models
        ecs_this_sample = set()
        models_not_in_parsing = []
        for mf in model_files:
            model_name = mf.replace('.xml', '')
            if model_name in model_to_ecs:
                ecs_this_sample.update(model_to_ecs[model_name])
            else:
                models_not_in_parsing.append(model_name)
        
        if samples_processed < 3:
            log(f"  {sample_id} ({xlsx_notation}): {len(model_files)} models, "
                f"{len(ecs_this_sample)} unique ECs"
                + (f", {len(models_not_in_parsing)} models not parsed" if models_not_in_parsing else ""))
        
        # Add rows to output
        for ec in ecs_this_sample:
            rows_out.append((sample_id, ec))
        
        samples_processed += 1
    
    log(f"\nSamples processed: {samples_processed}")
    log(f"Samples skipped: {len(samples_skipped)}")
    if samples_skipped:
        for s, x, r in samples_skipped[:10]:
            log(f"  {s} ({x}): {r}")
    
    # Write output
    log(f"\nWriting {len(rows_out)} (sample, EC) rows to {OUTPUT_TSV}")
    with open(OUTPUT_TSV, 'w') as f:
        f.write("sample_id\tec_number\n")
        for sample_id, ec in rows_out:
            f.write(f"{sample_id}\t{ec}\n")
    
    log("\n=== SUMMARY ===")
    log(f"Samples with AGORA EC sets: {samples_processed}")
    log(f"Total (sample, EC) pairs: {len(rows_out)}")
    log(f"Average ECs per sample: {len(rows_out) / samples_processed:.1f}")
    log(f"Output: {OUTPUT_TSV}")

if __name__ == "__main__":
    main()
