#!/usr/bin/env python3
"""
Step 5: Per-metabolite focused concordance analysis.

Markers updated based on Step 5b re-verification — corrected EC codes
that are present in BOTH HUMAnN3 and AGORA databases.

Outputs:
  concordance/metabolite_concordance.tsv
  concordance/ec_spearman_PRE_FBA.tsv  (pre-FBA Spearman, EC abundance vs AGORA presence)
  concordance/metabolite_spearman_PRE_FBA.tsv (aggregated)
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from collections import defaultdict

AGORA_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_agora_ecs_filtered.tsv"
HUMANN3_LONG_TSV = "/scratch/amans/concern3_analysis/per_sample/sample_humann3_ecs.tsv"
CRC_EC_CPM = "/scratch/amans/humann3_work/ec_grouped/crc_ec_cpm.tsv"
HEALTHY_EC_CPM = "/scratch/amans/humann3_work/ec_grouped/healthy_ec_cpm.tsv"
USABLE_SAMPLES = "/scratch/amans/concern3_analysis/per_sample/usable_samples.txt"
MAPPING_TSV = "/scratch/amans/bsim/sample_mapping_tab.tsv"
OUT_DIR = "/scratch/amans/concern3_analysis/concordance"

# VERIFIED MARKER ECs — present in BOTH HUMAnN3 and AGORA databases
METABOLITE_ECS = {
    'Butyrate':   ['2.7.2.7', '2.8.3.8', '2.8.3.9'],
    'Acetate':    ['2.7.2.1', '2.3.1.8'],
    'Propionate': ['5.4.99.2'],
    'Succinate':  ['1.3.5.2', '1.3.99.1'],
    'H2S':        ['4.4.1.1'],
    'TMA':        ['1.6.6.9'],
    'Putrescine': ['4.1.1.17', '4.1.1.19'],
}

def load_agora_per_sample():
    df = pd.read_csv(AGORA_TSV, sep='\t')
    result = defaultdict(set)
    for _, row in df.iterrows():
        result[row['sample_id']].add(row['ec_number'])
    return result

def load_humann3_per_sample_binary():
    df = pd.read_csv(HUMANN3_LONG_TSV, sep='\t')
    result = defaultdict(set)
    for _, row in df.iterrows():
        result[row['sample_id']].add(row['ec_number'])
    return result

def load_humann3_cpm_matrix(usable_samples):
    dfs = []
    for path in [CRC_EC_CPM, HEALTHY_EC_CPM]:
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df[~df.index.str.contains(r'\|', regex=True)]
        df.columns = [c.replace('_Abundance-CPM', '').replace('_Abundance-RPKs', '') for c in df.columns]
        df = df[~df.index.isin(['UNMAPPED', 'UNGROUPED'])]
        dfs.append(df)
    merged = pd.concat(dfs, axis=1)
    keep = [c for c in merged.columns if c in usable_samples]
    return merged[keep]

def main():
    print("=== Step 5c: Re-run with VERIFIED markers ===")
    
    with open(USABLE_SAMPLES) as f:
        usable = set(line.strip() for line in f if line.strip())
    
    agora_sets = load_agora_per_sample()
    humann3_sets = load_humann3_per_sample_binary()
    h_cpm = load_humann3_cpm_matrix(usable)
    
    df_map = pd.read_csv(MAPPING_TSV, sep='\t')
    sample_to_group = dict(zip(df_map['sample_id'], df_map['group']))
    
    common_samples = sorted(set(agora_sets) & set(humann3_sets) & set(h_cpm.columns))
    print(f"Samples in analysis: {len(common_samples)}")
    
    # === PART A: Per-metabolite concordance ===
    concordance_rows = []
    for metabolite, ecs in METABOLITE_ECS.items():
        ec_set = set(ecs)
        for group in ['CRC', 'Healthy']:
            samples_group = [s for s in common_samples if sample_to_group.get(s) == group]
            both = agora_only = humann3_only = neither = 0
            for s in samples_group:
                a_has = bool(agora_sets[s] & ec_set)
                h_has = bool(humann3_sets[s] & ec_set)
                if a_has and h_has: both += 1
                elif a_has: agora_only += 1
                elif h_has: humann3_only += 1
                else: neither += 1
            n = len(samples_group)
            concordance_rows.append({
                'metabolite': metabolite,
                'group': group,
                'n_samples': n,
                'both': both,
                'agora_only': agora_only,
                'humann3_only': humann3_only,
                'neither': neither,
                'pct_concordant': round(100 * both / n, 1) if n else 0,
                'pct_agora_supported_by_humann3': round(100 * both / (both + agora_only), 1) if (both + agora_only) else 0,
            })
    
    df_conc = pd.DataFrame(concordance_rows)
    df_conc.to_csv(f"{OUT_DIR}/metabolite_concordance.tsv", sep='\t', index=False)
    
    # === PART B: Spearman per EC (pre-FBA) ===
    spearman_ec_rows = []
    for metabolite, ecs in METABOLITE_ECS.items():
        for ec in ecs:
            if ec not in h_cpm.index:
                continue
            vec_h = []
            vec_a = []
            for s in common_samples:
                vec_h.append(h_cpm.at[ec, s])
                vec_a.append(1 if ec in agora_sets[s] else 0)
            vec_h = np.array(vec_h)
            vec_a = np.array(vec_a)
            if vec_h.std() == 0 or vec_a.std() == 0:
                rho, pval = np.nan, np.nan
            else:
                rho, pval = spearmanr(vec_h, vec_a)
            spearman_ec_rows.append({
                'metabolite': metabolite,
                'ec': ec,
                'n_samples': len(common_samples),
                'humann3_nonzero_n': int((vec_h > 0).sum()),
                'agora_present_n': int(vec_a.sum()),
                'spearman_rho': round(rho, 3) if not np.isnan(rho) else None,
                'pvalue': f"{pval:.2e}" if not np.isnan(pval) else None,
            })
    df_sp = pd.DataFrame(spearman_ec_rows)
    df_sp.to_csv(f"{OUT_DIR}/ec_spearman_PRE_FBA.tsv", sep='\t', index=False)
    
    spearman_met_rows = []
    for metabolite in METABOLITE_ECS:
        rhos = [r['spearman_rho'] for r in spearman_ec_rows if r['metabolite'] == metabolite and r['spearman_rho'] is not None]
        if rhos:
            mean_rho = np.mean(rhos)
            max_rho = max(rhos, key=abs)
        else:
            mean_rho = None
            max_rho = None
        spearman_met_rows.append({
            'metabolite': metabolite,
            'n_marker_ecs_used': len(rhos),
            'mean_rho': round(mean_rho, 3) if mean_rho is not None else None,
            'max_abs_rho': round(max_rho, 3) if max_rho is not None else None,
        })
    df_met = pd.DataFrame(spearman_met_rows)
    df_met.to_csv(f"{OUT_DIR}/metabolite_spearman_PRE_FBA.tsv", sep='\t', index=False)
    
    print("\n=== METABOLITE CONCORDANCE v2 ===")
    print(df_conc.to_string(index=False))
    print("\n=== EC SPEARMAN v2 (pre-FBA) ===")
    print(df_sp.to_string(index=False))
    print("\n=== METABOLITE SPEARMAN AGGREGATE v2 ===")
    print(df_met.to_string(index=False))

if __name__ == "__main__":
    main()
