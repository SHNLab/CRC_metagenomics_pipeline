#!/usr/bin/env python3
"""
Step 9: Build MAG summary tables for response letter.

Combines CheckM2 quality + GTDB-Tk taxonomy + bin-to-sample mapping
into analysis-ready tables.

Outputs:
  - data/02_checkm2/mag_summary.tsv (per-MAG combined info)
  - data/02_checkm2/mag_quality_summary.tsv (cohort-level stats)
  - data/03_gtdbtk/taxonomy_distribution.tsv (family/genus/species counts)
"""

import os
import pandas as pd
import numpy as np
import re

BASE = "/scratch/amans/concern3_analysis"
DELIV = f"{BASE}/deliverable"
CHECKM2_CRC = f"{DELIV}/data/02_checkm2/crc_checkm2_quality.tsv"
CHECKM2_HEALTHY = f"{DELIV}/data/02_checkm2/healthy_checkm2_quality.tsv"
GTDB_FILE = f"{DELIV}/data/03_gtdbtk/gtdbtk.bac120.summary.tsv"
KEEPERS_CRC = f"{DELIV}/data/02_checkm2/keepers_crc.txt"
KEEPERS_HEALTHY = f"{DELIV}/data/02_checkm2/keepers_healthy.txt"
OUT_MAG_SUMMARY = f"{DELIV}/data/02_checkm2/mag_summary.tsv"
OUT_QUALITY_SUMMARY = f"{DELIV}/data/02_checkm2/mag_quality_summary.tsv"
OUT_TAX_DIST = f"{DELIV}/data/03_gtdbtk/taxonomy_distribution.tsv"

def extract_tax_level(s, prefix):
    """Extract taxonomy level (e.g., 'f__' for family) from GTDB string."""
    if not s or pd.isna(s):
        return None
    for part in str(s).split(';'):
        part = part.strip()
        if part.startswith(prefix):
            name = part[len(prefix):].strip()
            return name if name else None
    return None

def extract_sample_id(bin_id):
    """Extract SRR ID from bin ID like SRR8865572.bin_1 or SRR8865572_bin_1."""
    m = re.match(r'(SRR\d+)', str(bin_id))
    return m.group(1) if m else None

def main():
    print("=== Step 9: MAG summary tables ===")
    
    # Load CheckM2
    crc_qc = pd.read_csv(CHECKM2_CRC, sep='\t')
    hlt_qc = pd.read_csv(CHECKM2_HEALTHY, sep='\t')
    crc_qc['cohort'] = 'CRC'
    hlt_qc['cohort'] = 'Healthy'
    qc = pd.concat([crc_qc, hlt_qc], ignore_index=True)
    print(f"Loaded {len(qc)} total bins from CheckM2")
    
    # Load keepers
    with open(KEEPERS_CRC) as f:
        keepers_crc = set(line.strip() for line in f if line.strip())
    with open(KEEPERS_HEALTHY) as f:
        keepers_hlt = set(line.strip() for line in f if line.strip())
    all_keepers = keepers_crc | keepers_hlt
    
    qc['is_keeper'] = qc['Name'].isin(all_keepers)
    print(f"Keepers: {qc['is_keeper'].sum()} ({len(keepers_crc)} CRC, {len(keepers_hlt)} Healthy)")
    
    # Load GTDB-Tk
    gtdb = pd.read_csv(GTDB_FILE, sep='\t')
    print(f"GTDB-Tk classified: {len(gtdb)} MAGs")
    
    # Split taxonomy into columns
    for level, prefix in [('domain', 'd__'), ('phylum', 'p__'), ('class', 'c__'),
                          ('order', 'o__'), ('family', 'f__'), ('genus', 'g__'), ('species', 's__')]:
        gtdb[level] = gtdb['classification'].apply(lambda s: extract_tax_level(s, prefix))
    
    # Merge
    merged = qc[qc['is_keeper']].merge(
        gtdb[['user_genome', 'classification', 'domain', 'phylum', 'class', 'order',
              'family', 'genus', 'species']],
        left_on='Name', right_on='user_genome', how='left'
    )
    
    # Add sample_id
    merged['sample_id'] = merged['Name'].apply(extract_sample_id)
    
    # Reorganize columns
    keep_cols = ['Name', 'sample_id', 'cohort', 'Completeness', 'Contamination',
                 'Completeness_Model_Used', 'Genome_Size', 'GC_Content', 'classification',
                 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    keep_cols = [c for c in keep_cols if c in merged.columns]
    mag_summary = merged[keep_cols].copy()
    mag_summary = mag_summary.rename(columns={'Name': 'bin_id'})
    mag_summary.to_csv(OUT_MAG_SUMMARY, sep='\t', index=False)
    print(f"\nSaved: {OUT_MAG_SUMMARY}")
    print(f"  {len(mag_summary)} keeper MAGs with combined quality + taxonomy")
    
    # Quality summary per cohort
    quality_rows = []
    for cohort in ['CRC', 'Healthy']:
        sub = qc[qc['cohort'] == cohort]
        sub_keepers = sub[sub['is_keeper']]
        quality_rows.append({
            'cohort': cohort,
            'total_bins': len(sub),
            'medium_quality_keepers': len(sub_keepers),
            'pct_passing': round(100 * len(sub_keepers) / len(sub), 1),
            'mean_completeness': round(sub_keepers['Completeness'].mean(), 1),
            'median_completeness': round(sub_keepers['Completeness'].median(), 1),
            'mean_contamination': round(sub_keepers['Contamination'].mean(), 2),
            'median_contamination': round(sub_keepers['Contamination'].median(), 2),
            'high_quality_count': int(((sub_keepers['Completeness'] >= 90) &
                                       (sub_keepers['Contamination'] < 5)).sum()),
        })
    
    df_quality = pd.DataFrame(quality_rows)
    df_quality.to_csv(OUT_QUALITY_SUMMARY, sep='\t', index=False)
    print(f"\nSaved: {OUT_QUALITY_SUMMARY}")
    print(df_quality.to_string(index=False))
    
    # Taxonomy distribution
    tax_rows = []
    for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        for cohort in ['CRC', 'Healthy']:
            sub = mag_summary[mag_summary['cohort'] == cohort]
            n_classified = sub[level].notna().sum()
            n_unique = sub[level].nunique()
            tax_rows.append({
                'taxonomic_level': level,
                'cohort': cohort,
                'n_mags_classified': int(n_classified),
                'n_unique_taxa': int(n_unique),
                'pct_resolution': round(100 * n_classified / len(sub), 1) if len(sub) else 0,
            })
    df_tax = pd.DataFrame(tax_rows)
    df_tax.to_csv(OUT_TAX_DIST, sep='\t', index=False)
    print(f"\nSaved: {OUT_TAX_DIST}")
    print(df_tax.to_string(index=False))
    
    # Headline summary
    print("\n=== HEADLINE ===")
    print(f"Total MAGs assembled: {len(qc)}")
    print(f"Quality-filtered keepers: {qc['is_keeper'].sum()} (>70% complete, <10% contam)")
    print(f"  - CRC: {len(keepers_crc)}")
    print(f"  - Healthy: {len(keepers_hlt)}")
    print(f"GTDB-Tk species-level resolution: {(mag_summary['species'].notna() & (mag_summary['species'] != '')).sum()}/{len(mag_summary)} ({100*(mag_summary['species'].notna() & (mag_summary['species'] != '')).sum()/len(mag_summary):.1f}%)")
    print(f"Patients with at least 1 keeper MAG: {mag_summary['sample_id'].nunique()}")

if __name__ == "__main__":
    main()
