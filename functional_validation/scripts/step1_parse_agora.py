#!/usr/bin/env python3
"""
Step 1: Parse AGORA .xml models and extract EC numbers per model.

Output: /scratch/amans/concern3_analysis/parsing/agora_model_ecs.tsv
Columns: model_name, ec_number, n_reactions_with_this_ec

Method: Uses cobrapy to read SBML models, extracts EC annotations from
reaction's annotation dict (key 'ec-code'), splits multi-EC entries
(e.g. "2.7.1.167 2.7.7.70") on whitespace.
"""

import os
import sys
import glob
from collections import defaultdict
import cobra
from cobra.io import read_sbml_model

MODELS_DIR = "/scratch/amans/bsim/models_fixed"
OUTPUT_TSV = "/scratch/amans/concern3_analysis/parsing/agora_model_ecs.tsv"
LOG_FILE = "/scratch/amans/concern3_analysis/logs/step1_parse_agora.log"

def log(msg, also_print=True):
    with open(LOG_FILE, 'a') as f:
        f.write(msg + '\n')
    if also_print:
        print(msg)

def extract_ec_from_reaction(reaction):
    """Extract EC numbers — split multi-EC entries on whitespace."""
    annot = reaction.annotation
    ec_codes = []
    
    if 'ec-code' in annot:
        val = annot['ec-code']
        if isinstance(val, list):
            for v in val:
                ec_codes.extend(str(v).split())
        else:
            ec_codes.extend(str(val).split())
    
    return [ec.strip() for ec in ec_codes if ec.strip()]

def parse_model(model_path):
    """Parse one AGORA SBML model, return dict of {ec: count_of_reactions}."""
    model_name = os.path.basename(model_path).replace('.xml', '')
    ec_counts = defaultdict(int)
    
    try:
        model = read_sbml_model(model_path)
    except Exception as e:
        log(f"  ERROR loading {model_name}: {e}")
        return model_name, None
    
    for rxn in model.reactions:
        ecs = extract_ec_from_reaction(rxn)
        for ec in ecs:
            ec_counts[ec] += 1
    
    return model_name, dict(ec_counts)

def main():
    open(LOG_FILE, 'w').close()
    log(f"=== Step 1: AGORA EC extraction ===")
    log(f"Models dir: {MODELS_DIR}")
    
    model_files = sorted(glob.glob(os.path.join(MODELS_DIR, "*.xml")))
    log(f"Found {len(model_files)} .xml models")
    
    if len(model_files) == 0:
        log("FATAL: no .xml files found")
        sys.exit(1)
    
    all_rows = []
    failed_models = []
    models_with_no_ecs = []
    
    for i, mp in enumerate(model_files, 1):
        if i % 20 == 0 or i == 1 or i == len(model_files):
            log(f"  Processing {i}/{len(model_files)}: {os.path.basename(mp)}")
        
        model_name, ec_counts = parse_model(mp)
        
        if ec_counts is None:
            failed_models.append(model_name)
            continue
        
        if len(ec_counts) == 0:
            models_with_no_ecs.append(model_name)
            continue
        
        for ec, count in ec_counts.items():
            all_rows.append((model_name, ec, count))
    
    log(f"\nWriting {len(all_rows)} (model, EC) rows to {OUTPUT_TSV}")
    with open(OUTPUT_TSV, 'w') as f:
        f.write("model_name\tec_number\tn_reactions_with_this_ec\n")
        for model_name, ec, count in all_rows:
            f.write(f"{model_name}\t{ec}\t{count}\n")
    
    unique_ecs = set(row[1] for row in all_rows)
    log(f"\n=== SUMMARY ===")
    log(f"Models processed successfully: {len(model_files) - len(failed_models) - len(models_with_no_ecs)}")
    log(f"Models with no EC annotations: {len(models_with_no_ecs)}")
    log(f"Failed to load: {len(failed_models)}")
    log(f"Unique ECs across all models: {len(unique_ecs)}")
    log(f"Total (model, EC) pairs: {len(all_rows)}")
    
    if failed_models:
        log(f"\nFailed models:")
        for fm in failed_models[:20]:
            log(f"  {fm}")
        if len(failed_models) > 20:
            log(f"  ... and {len(failed_models) - 20} more")

if __name__ == "__main__":
    main()
