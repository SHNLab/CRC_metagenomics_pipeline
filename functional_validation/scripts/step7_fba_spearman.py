#!/usr/bin/env python3
"""
Step 7: Post-FBA Spearman correlation between HUMAnN3 abundance and AGORA flux.

Replaces step7 with corrected EC list.
"""

import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

FBA_BASE = "/scratch/amans/bsim/FBA"
DIETS = ['FIBER', 'GLUTENFREE', 'HFD', 'PROTEIN', 'VEGAN']
USABLE_SAMPLES = "/scratch/amans/concern3_analysis/per_sample/usable_samples.txt"
MAPPING_TSV = "/scratch/amans/bsim/sample_mapping_tab.tsv"
CRC_EC_CPM = "/scratch/amans/humann3_work/ec_grouped/crc_ec_cpm.tsv"
HEALTHY_EC_CPM = "/scratch/amans/humann3_work/ec_grouped/healthy_ec_cpm.tsv"
CRC_PWY_CPM = "/scratch/amans/humann3_work/normalized/crc_pathabundance_cpm.tsv"
HEALTHY_PWY_CPM = "/scratch/amans/humann3_work/normalized/healthy_pathabundance_cpm.tsv"
OUT_DIR = "/scratch/amans/concern3_analysis/concordance"

# VERIFIED markers
METABOLITE_DEFS = {
    'Butyrate':   {'ecs': ['2.7.2.7', '2.8.3.8', '2.8.3.9'], 'pwys': ['PWY-5676', 'PWY-5677', 'CENTFERM-PWY', 'P163-PWY']},
    'Acetate':    {'ecs': ['2.7.2.1', '2.3.1.8'],            'pwys': ['P41-PWY', 'PWY-5100', 'PWY-7254']},
    'Propionate': {'ecs': ['5.4.99.2'],                       'pwys': ['P108-PWY', 'PWY-8086']},
    'Succinate':  {'ecs': ['1.3.5.2', '1.3.99.1'],            'pwys': ['PWY-5677']},
    'Hydrogen sulfide': {'ecs': ['4.4.1.1'],                  'pwys': []},
    'Trimethylamine':   {'ecs': ['1.6.6.9'],                  'pwys': []},
    'Putrescine': {'ecs': ['4.1.1.17', '4.1.1.19'],           'pwys': ['PWY-6305', 'PWY0-1221']},
}

def parse_scfa_csv(path):
    df = pd.read_csv(path)
    if len(df) == 0:
        return {}
    row = df.iloc[0]
    result = {}
    for i in range(1, 14):
        met_col = f"Metabolite_{i}"
        flux_col = f"WeightedFlux_{i}"
        if met_col in df.columns and flux_col in df.columns:
            met = str(row[met_col]).strip()
            try:
                flux = float(row[flux_col])
            except:
                flux = np.nan
            result[met] = flux
    return result

def patient_from_filename(filepath):
    fn = os.path.basename(filepath).replace('.csv', '')
    for sep in ['_di_p', '_hi_p']:
        idx = fn.rfind(sep)
        if idx >= 0:
            return sep[1:] + fn[idx + len(sep):]
    return None

def main():
    print("=== Step 7: Post-FBA Spearman correlation between HUMAnN3 abundance and AGORA flux.===")
    
    with open(USABLE_SAMPLES) as f:
        usable = set(line.strip() for line in f if line.strip())
    
    df_map = pd.read_csv(MAPPING_TSV, sep='\t')
    sample_to_xlsx = dict(zip(df_map['sample_id'], df_map['xlxs_notation']))
    xlsx_to_sample = {v: k for k, v in sample_to_xlsx.items()}
    
    # Load HUMAnN3 EC matrix
    dfs = []
    for path in [CRC_EC_CPM, HEALTHY_EC_CPM]:
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df[~df.index.str.contains(r'\|', regex=True)]
        df = df[~df.index.isin(['UNMAPPED', 'UNGROUPED'])]
        df.columns = [c.replace('_Abundance-CPM', '').replace('_Abundance-RPKs', '') for c in df.columns]
        dfs.append(df)
    h3_ec = pd.concat(dfs, axis=1)
    h3_ec = h3_ec[[c for c in h3_ec.columns if c in usable]]
    
    # Load HUMAnN3 pathway matrix
    dfs = []
    for path in [CRC_PWY_CPM, HEALTHY_PWY_CPM]:
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df[~df.index.str.contains(r'\|', regex=True)]
        df.columns = [c.replace('_Abundance-CPM', '').replace('_Abundance', '') for c in df.columns]
        dfs.append(df)
    h3_pwy = pd.concat(dfs, axis=1)
    h3_pwy = h3_pwy[[c for c in h3_pwy.columns if c in usable]]
    h3_pwy.index = h3_pwy.index.map(lambda x: x.split(':')[0].strip())
    
    all_results = []
    for diet in DIETS:
        scfa_dir = os.path.join(FBA_BASE, f"{diet}_FLUX")
        scfa_files = glob.glob(os.path.join(scfa_dir, "**", "SCFA*.csv"), recursive=True)
        
        diet_flux = {}
        for scfa in scfa_files:
            pid = patient_from_filename(scfa)
            if pid is None:
                continue
            sid = xlsx_to_sample.get(pid)
            if sid is None or sid not in usable:
                continue
            diet_flux[sid] = parse_scfa_csv(scfa)
        
        for met, defs in METABOLITE_DEFS.items():
            samples = sorted(set(diet_flux.keys()) & set(h3_ec.columns))
            if len(samples) < 5:
                continue
            
            vec_a = np.array([diet_flux[s].get(met, np.nan) for s in samples])
            
            ec_rows = [e for e in defs['ecs'] if e in h3_ec.index]
            vec_ec = h3_ec.loc[ec_rows, samples].sum(axis=0).values if ec_rows else np.full(len(samples), np.nan)
            
            pwy_rows = [p for p in defs['pwys'] if p in h3_pwy.index]
            vec_pwy = h3_pwy.loc[pwy_rows, samples].sum(axis=0).values if pwy_rows else np.full(len(samples), np.nan)
            
            mask_ec = ~(np.isnan(vec_a) | np.isnan(vec_ec))
            if mask_ec.sum() >= 5 and vec_a[mask_ec].std() > 0 and vec_ec[mask_ec].std() > 0:
                rho_ec, p_ec = spearmanr(vec_a[mask_ec], vec_ec[mask_ec])
            else:
                rho_ec, p_ec = np.nan, np.nan
            
            mask_pwy = ~(np.isnan(vec_a) | np.isnan(vec_pwy))
            if mask_pwy.sum() >= 5 and vec_a[mask_pwy].std() > 0 and vec_pwy[mask_pwy].std() > 0:
                rho_pwy, p_pwy = spearmanr(vec_a[mask_pwy], vec_pwy[mask_pwy])
            else:
                rho_pwy, p_pwy = np.nan, np.nan
            
            all_results.append({
                'diet': diet,
                'metabolite': met,
                'n_samples': len(samples),
                'n_ec_markers': len(ec_rows),
                'n_pwy_markers': len(pwy_rows),
                'rho_ec_vs_flux': round(rho_ec, 3) if not np.isnan(rho_ec) else None,
                'p_ec_vs_flux': f"{p_ec:.2e}" if not np.isnan(p_ec) else None,
                'rho_pwy_vs_flux': round(rho_pwy, 3) if not np.isnan(rho_pwy) else None,
                'p_pwy_vs_flux': f"{p_pwy:.2e}" if not np.isnan(p_pwy) else None,
            })
    
    df_r = pd.DataFrame(all_results)
    df_r.to_csv(f"{OUT_DIR}/spearman_per_diet_per_metabolite.tsv", sep='\t', index=False)
    
    # Aggregate per metabolite
    summary = []
    for met in METABOLITE_DEFS:
        sub = df_r[df_r['metabolite'] == met]
        rhos_ec = sub['rho_ec_vs_flux'].dropna().astype(float).tolist()
        rhos_pwy = sub['rho_pwy_vs_flux'].dropna().astype(float).tolist()
        summary.append({
            'metabolite': met,
            'n_diets_ec': len(rhos_ec),
            'median_rho_ec': round(np.median(rhos_ec), 3) if rhos_ec else None,
            'max_abs_rho_ec': round(max(rhos_ec, key=abs), 3) if rhos_ec else None,
            'n_diets_pwy': len(rhos_pwy),
            'median_rho_pwy': round(np.median(rhos_pwy), 3) if rhos_pwy else None,
            'max_abs_rho_pwy': round(max(rhos_pwy, key=abs), 3) if rhos_pwy else None,
        })
    df_s = pd.DataFrame(summary)
    df_s.to_csv(f"{OUT_DIR}/spearman_summary_per_metabolite.tsv", sep='\t', index=False)
    
    print("\n=== FULL TABLE v2 ===")
    print(df_r.to_string(index=False))
    print("\n=== SUMMARY v2 ===")
    print(df_s.to_string(index=False))

if __name__ == "__main__":
    main()
