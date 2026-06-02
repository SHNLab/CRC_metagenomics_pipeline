#!/usr/bin/env python3
"""
Step 3: Per-sample HUMAnN3 EC presence.

For each usable sample (129):
  - Read its column from crc_ec_cpm.tsv or healthy_ec_cpm.tsv
  - Keep only unstratified EC rows (no '|' in row name)
  - An EC is "present" if abundance > 0 CPM
  - Output: long-format (sample_id, ec_number)

Output: /scratch/amans/concern3_analysis/per_sample/sample_humann3_ecs.tsv
"""

import os
import sys
import pandas as pd

CRC_EC_TSV = "/scratch/amans/humann3_work/ec_grouped/crc_ec_cpm.tsv"
HEALTHY_EC_TSV = "/scratch/amans/humann3_work/ec_grouped/healthy_ec_cpm.tsv"
USABLE_SAMPLES = "/scratch/amans/concern3_analysis/per_sample/usable_samples.txt"
OUTPUT_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_humann3_ecs.tsv"
LOG_FILE = "/scratch/amans/concern3_analysis/logs/step3_build_sample_humann3_ec.log"

def log(msg, also_print=True):
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')
    if also_print:
        print(msg)

def process_cohort(ec_tsv, usable_samples):
    """Read EC matrix, return long-format DataFrame of (sample, ec) for unstratified ECs > 0."""
    log(f"Loading {ec_tsv}")
    df = pd.read_csv(ec_tsv, sep='\t', index_col=0)
    
    # Filter to unstratified rows (no '|' in index)
    unstratified = df[~df.index.str.contains(r'\|', regex=True)]
    log(f"  Total rows: {len(df)}, unstratified: {len(unstratified)}")
    
    # Filter to samples in usable list
    # Sample columns have format "SRR<id>_Abundance-CPM"
    rename_dict = {}
    for col in unstratified.columns:
        sample_id = col.replace('_Abundance-CPM', '').replace('_Abundance-RPKs', '')
        rename_dict[col] = sample_id
    unstratified = unstratified.rename(columns=rename_dict)
    
    keep_cols = [c for c in unstratified.columns if c in usable_samples]
    log(f"  Samples in cohort: {len(unstratified.columns)}, intersect with usable: {len(keep_cols)}")
    
    unstratified = unstratified[keep_cols]
    
    # Drop UNMAPPED and UNGROUPED rows
    unstratified = unstratified[~unstratified.index.isin(['UNMAPPED', 'UNGROUPED'])]
    log(f"  Rows after dropping UNMAPPED/UNGROUPED: {len(unstratified)}")
    
    # Convert to long format, keep only EC > 0 (presence)
    rows = []
    for sample in unstratified.columns:
        present_ecs = unstratified.index[unstratified[sample] > 0].tolist()
        for ec in present_ecs:
            rows.append((sample, ec))
    
    return rows

def main():
    open(LOG_FILE, 'w').close()
    log("=== Step 3: per-sample HUMAnN3 EC presence ===")
    
    # Load usable sample whitelist
    with open(USABLE_SAMPLES) as f:
        usable = set(line.strip() for line in f if line.strip())
    log(f"Usable samples: {len(usable)}")
    
    # Process both cohorts
    all_rows = []
    for label, path in [('CRC', CRC_EC_TSV), ('Healthy', HEALTHY_EC_TSV)]:
        log(f"\n--- {label} ---")
        rows = process_cohort(path, usable)
        log(f"  {label}: {len(rows)} (sample, EC) presence rows")
        all_rows.extend(rows)
    
    # Write
    log(f"\nWriting {len(all_rows)} rows to {OUTPUT_TSV}")
    with open(OUTPUT_TSV, 'w') as f:
        f.write("sample_id\tec_number\n")
        for sample, ec in all_rows:
            f.write(f"{sample}\t{ec}\n")
    
    # Summary
    samples_found = set(r[0] for r in all_rows)
    log("\n=== SUMMARY ===")
    log(f"Samples with HUMAnN3 EC data: {len(samples_found)}")
    log(f"Total (sample, EC) pairs: {len(all_rows)}")
    log(f"Average ECs per sample: {len(all_rows)/len(samples_found):.1f}")
    log(f"Output: {OUTPUT_TSV}")

if __name__ == "__main__":
    main()
