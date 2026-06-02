#!/usr/bin/env python3
"""
Step 5b: Re-verify marker ECs (FAST version).

Loads each AGORA model ONCE, searches all keywords inside that model,
accumulates results across models. ~10-15 min total.

Outputs: parsing/marker_ec_verification.tsv
"""

import os
import glob
import pandas as pd
import cobra
from cobra.io import read_sbml_model

AGORA_DIR = "/scratch/amans/bsim/models_fixed"
CRC_EC = "/scratch/amans/humann3_work/ec_grouped/crc_ec_cpm.tsv"
HEALTHY_EC = "/scratch/amans/humann3_work/ec_grouped/healthy_ec_cpm.tsv"
OUT_TSV = "/scratch/amans/concern3_analysis/parsing/marker_ec_verification.tsv"

SEARCHES = [
    ('Butyrate', 'butyrate kinase', ['BUTYRATE KINASE', 'BUTK'], ['2.7.2.7']),
    ('Butyrate', 'butyryl-CoA:acetate CoA-transferase', ['BUTYRYL-COA', 'ACETATE COA-TRANSFERASE'], ['2.8.3.8', '2.8.3.9']),
    ('Acetate', 'acetate kinase', ['ACETATE KINASE', 'ACKR'], ['2.7.2.1']),
    ('Acetate', 'phosphate acetyltransferase', ['PHOSPHATE ACETYLTRANSFERASE', 'PTAR'], ['2.3.1.8']),
    ('Propionate', 'methylmalonyl-CoA mutase', ['METHYLMALONYL-COA MUTASE'], ['5.4.99.2']),
    ('Propionate', 'methylmalonyl-CoA decarboxylase', ['METHYLMALONYL-COA DECARBOXYLASE'], ['4.1.1.41', '7.2.4.3']),
    ('Succinate', 'succinate dehydrogenase', ['SUCCINATE DEHYDROGENASE', 'SUCD'], ['1.3.5.1', '1.3.5.2', '1.3.99.1']),
    ('Succinate', 'fumarate reductase', ['FUMARATE REDUCTASE', 'FRD'], ['1.3.5.4', '1.3.5.1', '1.3.99.1']),
    ('H2S', 'dissimilatory sulfite reductase', ['DISSIMILATORY SULFITE REDUCTASE', 'DSRA'], ['1.8.99.5', '1.8.99.3', '1.8.99.1']),
    ('H2S', 'cysteine desulfhydrase', ['CYSTEINE DESULFHYDRASE'], ['4.4.1.1', '4.4.1.28']),
    ('TMA', 'choline TMA-lyase', ['CHOLINE TMA-LYASE', 'CUTC'], ['4.3.99.4']),
    ('TMA', 'TMAO reductase', ['TMAO REDUCTASE', 'TRIMETHYLAMINE'], ['1.7.2.3']),
    ('Putrescine', 'ornithine decarboxylase', ['ORNITHINE DECARBOXYLASE'], ['4.1.1.17']),
    ('Putrescine', 'arginine decarboxylase', ['ARGININE DECARBOXYLASE'], ['4.1.1.19']),
]

def load_humann3_ecs():
    ecs = set()
    for path in [CRC_EC, HEALTHY_EC]:
        df = pd.read_csv(path, sep='\t', index_col=0, usecols=[0])
        for idx in df.index:
            if '|' not in str(idx) and str(idx) not in ['UNMAPPED', 'UNGROUPED']:
                ecs.add(str(idx).strip())
    return ecs

def main():
    print("Loading HUMAnN3 EC list...")
    humann3_ecs = load_humann3_ecs()
    print(f"  {len(humann3_ecs)} unique ECs in HUMAnN3")
    
    agora_ecs_per_search = {i: set() for i in range(len(SEARCHES))}
    n_models_per_search = {i: 0 for i in range(len(SEARCHES))}
    
    model_files = sorted(glob.glob(os.path.join(AGORA_DIR, "*.xml")))
    print(f"Scanning {len(model_files)} models...")
    
    for n_done, mf in enumerate(model_files, 1):
        try:
            m = read_sbml_model(mf)
        except Exception as e:
            print(f"  Skip {os.path.basename(mf)}: {e}")
            continue
        
        rxn_info = [(r.id.upper(), (r.name or '').upper(), r.annotation.get('ec-code'))
                    for r in m.reactions]
        
        for i, (metab, enzyme, keywords, _lit_ecs) in enumerate(SEARCHES):
            matched_this_model = False
            for rid, rname, ec in rxn_info:
                for kw in keywords:
                    if kw in rid or kw in rname:
                        if ec:
                            if isinstance(ec, list):
                                for e in ec:
                                    agora_ecs_per_search[i].update(str(e).split())
                            else:
                                agora_ecs_per_search[i].update(str(ec).split())
                        matched_this_model = True
                        break
            if matched_this_model:
                n_models_per_search[i] += 1
        
        if n_done % 20 == 0 or n_done == len(model_files):
            print(f"  Processed {n_done}/{len(model_files)} models")
    
    rows = []
    for i, (metab, enzyme, keywords, lit_ecs) in enumerate(SEARCHES):
        agora_ecs = agora_ecs_per_search[i]
        n_models = n_models_per_search[i]
        
        for lec in lit_ecs:
            rows.append({
                'metabolite': metab,
                'enzyme': enzyme,
                'ec_source': 'literature',
                'ec': lec,
                'in_humann3': lec in humann3_ecs,
                'in_agora_via_ec_annotation': lec in agora_ecs,
                'n_agora_models_with_enzyme': n_models,
            })
        
        extra = agora_ecs - set(lit_ecs)
        for e in sorted(extra):
            rows.append({
                'metabolite': metab,
                'enzyme': enzyme,
                'ec_source': 'agora_by_name',
                'ec': e,
                'in_humann3': e in humann3_ecs,
                'in_agora_via_ec_annotation': True,
                'n_agora_models_with_enzyme': n_models,
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep='\t', index=False)
    
    print(f"\n=== SUMMARY ===")
    print(df.to_string(index=False))
    print(f"\nSaved: {OUT_TSV}")

if __name__ == "__main__":
    main()
