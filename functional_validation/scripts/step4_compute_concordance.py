#!/usr/bin/env python3
"""
Step 4: Per-sample EC concordance between AGORA and HUMAnN3.

Inputs:
  - per_sample/sample_agora_ecs_filtered.tsv (long format)
  - per_sample/sample_humann3_ecs.tsv (long format)
  - bsim/sample_mapping_tab.tsv (for cohort assignment)

Output:
  - concordance/concordance_per_sample.tsv
    Columns: sample_id, group, n_agora, n_humann3, n_intersection,
             n_agora_only, n_humann3_only, jaccard, pct_agora_explained, pct_humann3_explained
"""

import pandas as pd
from collections import defaultdict

AGORA_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_agora_ecs_filtered.tsv"
HUMANN3_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_humann3_ecs.tsv"
MAPPING_TSV = "/scratch/amans/bsim/sample_mapping_tab.tsv"
OUTPUT_TSV = "/scratch/amans/concern3_analysis/concordance/concordance_per_sample.tsv"
LOG_FILE = "/scratch/amans/concern3_analysis/logs/step4_compute_concordance.log"

def log(msg, also_print=True):
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')
    if also_print:
        print(msg)

def load_per_sample(tsv_path):
    """Load long-format (sample, EC) and return dict: sample -> set of ECs."""
    df = pd.read_csv(tsv_path, sep='\t')
    result = defaultdict(set)
    for _, row in df.iterrows():
        result[row['sample_id']].add(row['ec_number'])
    return result

def main():
    open(LOG_FILE, 'w').close()
    log("=== Step 4: Per-sample concordance ===")
    
    log("Loading AGORA per-sample ECs...")
    agora_ecs = load_per_sample(AGORA_TSV)
    log(f"  {len(agora_ecs)} samples")
    
    log("Loading HUMAnN3 per-sample ECs...")
    humann3_ecs = load_per_sample(HUMANN3_TSV)
    log(f"  {len(humann3_ecs)} samples")
    
    log("Loading sample mapping for cohort labels...")
    df_map = pd.read_csv(MAPPING_TSV, sep='\t')
    sample_to_group = dict(zip(df_map['sample_id'], df_map['group']))
    
    # Common samples
    common = set(agora_ecs.keys()) & set(humann3_ecs.keys())
    log(f"  Samples with both AGORA and HUMAnN3 data: {len(common)}")
    
    # Compute concordance per sample
    results = []
    for sample in sorted(common):
        a = agora_ecs[sample]
        h = humann3_ecs[sample]
        inter = a & h
        only_a = a - h
        only_h = h - a
        union = a | h
        
        jaccard = len(inter) / len(union) if union else 0.0
        pct_agora = len(inter) / len(a) * 100 if a else 0.0
        pct_humann3 = len(inter) / len(h) * 100 if h else 0.0
        
        results.append({
            'sample_id': sample,
            'group': sample_to_group.get(sample, 'UNKNOWN'),
            'n_agora': len(a),
            'n_humann3': len(h),
            'n_intersection': len(inter),
            'n_agora_only': len(only_a),
            'n_humann3_only': len(only_h),
            'jaccard': jaccard,
            'pct_agora_explained': pct_agora,
            'pct_humann3_explained': pct_humann3,
        })
    
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_TSV, sep='\t', index=False)
    
    log(f"\nOutput: {OUTPUT_TSV}")
    log("\n=== SUMMARY by cohort ===")
    summary = df.groupby('group').agg(
        n_samples=('sample_id', 'count'),
        jaccard_mean=('jaccard', 'mean'),
        jaccard_median=('jaccard', 'median'),
        pct_agora_mean=('pct_agora_explained', 'mean'),
        pct_humann3_mean=('pct_humann3_explained', 'mean'),
    ).round(3)
    log(str(summary))
    
    log("\n=== OVERALL ===")
    log(f"Mean Jaccard: {df['jaccard'].mean():.3f}")
    log(f"Median Jaccard: {df['jaccard'].median():.3f}")
    log(f"Mean % AGORA ECs explained by HUMAnN3: {df['pct_agora_explained'].mean():.1f}%")
    log(f"Mean % HUMAnN3 ECs explained by AGORA: {df['pct_humann3_explained'].mean():.1f}%")

if __name__ == "__main__":
    main()
