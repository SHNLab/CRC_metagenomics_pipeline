#!/usr/bin/env python3
"""
Generate Concern 3 figures locally on MacBook.
Reads pre-computed tables from data/ folder.
Produces 5 figures as PNG (300 dpi) and PDF in figures/ folder.

Usage:
    cd ~/Documents/deliverable
    python3 scripts/plot_concern3_figures.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DELIVERABLE = os.path.dirname(SCRIPT_DIR)
DATA = os.path.join(DELIVERABLE, "data")
OUT = os.path.join(DELIVERABLE, "figures")
os.makedirs(OUT, exist_ok=True)

sns.set_style("whitegrid")
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['axes.labelsize'] = 11
SAVE_DPI = 300

CRC_COLOR = "#c0392b"
HEALTHY_COLOR = "#2980b9"

def fig1_jaccard():
    df = pd.read_csv(os.path.join(DATA, "05_concordance", "concordance_per_sample.tsv"), sep='\t')
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    sns.violinplot(data=df, x='group', y='jaccard', ax=axes[0],
                   palette=[CRC_COLOR, HEALTHY_COLOR], inner='quartile')
    axes[0].set_ylabel("Jaccard index"); axes[0].set_xlabel("")
    axes[0].set_title("Per-sample Jaccard concordance"); axes[0].set_ylim(0, 0.6)
    sns.violinplot(data=df, x='group', y='pct_agora_explained', ax=axes[1],
                   palette=[CRC_COLOR, HEALTHY_COLOR], inner='quartile')
    axes[1].set_ylabel("% AGORA ECs supported by HUMAnN3"); axes[1].set_xlabel("")
    axes[1].set_title("AGORA-side support"); axes[1].set_ylim(0, 100)
    sns.violinplot(data=df, x='group', y='pct_humann3_explained', ax=axes[2],
                   palette=[CRC_COLOR, HEALTHY_COLOR], inner='quartile')
    axes[2].set_ylabel("% HUMAnN3 ECs covered by AGORA"); axes[2].set_xlabel("")
    axes[2].set_title("HUMAnN3-side coverage"); axes[2].set_ylim(0, 100)
    fig.suptitle("Figure 1. Pre-FBA structural concordance (n=129)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig1_jaccard_concordance.png"), dpi=SAVE_DPI, bbox_inches='tight')
    plt.savefig(os.path.join(OUT, "fig1_jaccard_concordance.pdf"), bbox_inches='tight')
    plt.close()
    print("Saved fig1_jaccard_concordance.{png,pdf}")

def fig2_metabolite_concordance():
    df = pd.read_csv(os.path.join(DATA, "05_concordance", "metabolite_concordance.tsv"), sep='\t')
    met_order = ['Acetate', 'Butyrate', 'Propionate', 'Putrescine', 'Succinate', 'H2S', 'TMA']
    df['metabolite'] = pd.Categorical(df['metabolite'], categories=met_order, ordered=True)
    df = df.sort_values('metabolite')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, group in zip(axes, ['CRC', 'Healthy']):
        sub = df[df['group'] == group].copy()
        n = sub['n_samples'].iloc[0]
        sub['pct_both'] = 100 * sub['both'] / n
        sub['pct_agora_only'] = 100 * sub['agora_only'] / n
        sub['pct_humann3_only'] = 100 * sub['humann3_only'] / n
        sub['pct_neither'] = 100 * sub['neither'] / n
        x = np.arange(len(sub))
        ax.bar(x, sub['pct_both'], label='Both', color='#27ae60')
        ax.bar(x, sub['pct_agora_only'], bottom=sub['pct_both'], label='AGORA only', color='#e67e22')
        ax.bar(x, sub['pct_humann3_only'],
               bottom=sub['pct_both'] + sub['pct_agora_only'], label='HUMAnN3 only', color='#8e44ad')
        ax.bar(x, sub['pct_neither'],
               bottom=sub['pct_both'] + sub['pct_agora_only'] + sub['pct_humann3_only'],
               label='Neither', color='#7f8c8d')
        ax.set_xticks(x); ax.set_xticklabels(sub['metabolite'], rotation=45, ha='right')
        ax.set_ylabel("% of samples"); ax.set_title(f"{group} (n={n})"); ax.set_ylim(0, 105)
        if group == 'Healthy':
            ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1))
    fig.suptitle("Figure 2. Per-metabolite marker EC concordance", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig2_metabolite_concordance.png"), dpi=SAVE_DPI, bbox_inches='tight')
    plt.savefig(os.path.join(OUT, "fig2_metabolite_concordance.pdf"), bbox_inches='tight')
    plt.close()
    print("Saved fig2_metabolite_concordance.{png,pdf}")

def fig3_spearman_heatmap():
    df = pd.read_csv(os.path.join(DATA, "06_spearman", "spearman_per_diet_per_metabolite.tsv"), sep='\t')
    diet_order = ['FIBER', 'GLUTENFREE', 'HFD', 'PROTEIN', 'VEGAN']
    met_order = ['Acetate', 'Butyrate', 'Propionate', 'Putrescine', 'Succinate', 'Hydrogen sulfide', 'Trimethylamine']
    df_ec = df.pivot(index='metabolite', columns='diet', values='rho_ec_vs_flux')
    df_pwy = df.pivot(index='metabolite', columns='diet', values='rho_pwy_vs_flux')
    df_ec = df_ec.reindex(index=met_order, columns=diet_order)
    df_pwy = df_pwy.reindex(index=met_order, columns=diet_order)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    sns.heatmap(df_ec, annot=True, fmt='.2f', cmap=cmap, center=0, vmin=-0.5, vmax=0.5,
                cbar_kws={'label': 'Spearman rho'}, ax=axes[0], linewidths=0.5)
    axes[0].set_title("EC-sum vs AGORA flux"); axes[0].set_xlabel("Diet"); axes[0].set_ylabel("Metabolite")
    sns.heatmap(df_pwy, annot=True, fmt='.2f', cmap=cmap, center=0, vmin=-0.5, vmax=0.5,
                cbar_kws={'label': 'Spearman rho'}, ax=axes[1], linewidths=0.5)
    axes[1].set_title("Pathway-sum vs AGORA flux"); axes[1].set_xlabel("Diet"); axes[1].set_ylabel("")
    fig.suptitle("Figure 3. Diet-dependent Spearman correlation (n=129)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig3_spearman_heatmap.png"), dpi=SAVE_DPI, bbox_inches='tight')
    plt.savefig(os.path.join(OUT, "fig3_spearman_heatmap.pdf"), bbox_inches='tight')
    plt.close()
    print("Saved fig3_spearman_heatmap.{png,pdf}")

def fig4_checkm2_quality():
    crc = pd.read_csv(os.path.join(DATA, "02_checkm2", "crc_checkm2_quality.tsv"), sep='\t')
    hlt = pd.read_csv(os.path.join(DATA, "02_checkm2", "healthy_checkm2_quality.tsv"), sep='\t')
    crc['cohort'] = 'CRC'; hlt['cohort'] = 'Healthy'
    df = pd.concat([crc, hlt], ignore_index=True)
    fig, ax = plt.subplots(figsize=(9, 7))
    for cohort, color in [('CRC', CRC_COLOR), ('Healthy', HEALTHY_COLOR)]:
        sub = df[df['cohort'] == cohort]
        ax.scatter(sub['Completeness'], sub['Contamination'], alpha=0.6, s=30,
                   color=color, label=f"{cohort} (n={len(sub)})", edgecolor='white', linewidth=0.5)
    ax.axhline(10, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(70, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(90, color='black', linestyle=':', linewidth=1, alpha=0.5)
    ax.axhline(5, color='black', linestyle=':', linewidth=1, alpha=0.5)
    ax.text(95, 2.5, "High\nquality", ha='center', va='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
    ax.text(82, 7.5, "Medium\nquality", ha='center', va='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    ax.set_xlabel("Completeness (%)"); ax.set_ylabel("Contamination (%)")
    ax.set_xlim(0, 105); ax.set_ylim(0, min(50, df['Contamination'].max() * 1.05))
    ax.legend(loc='upper left')
    ax.set_title("Figure 4. CheckM2 quality of MAGs\n(746 total bins, 226 keepers)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig4_checkm2_quality.png"), dpi=SAVE_DPI, bbox_inches='tight')
    plt.savefig(os.path.join(OUT, "fig4_checkm2_quality.pdf"), bbox_inches='tight')
    plt.close()
    print("Saved fig4_checkm2_quality.{png,pdf}")

def fig5_gtdbtk_taxonomy():
    df = pd.read_csv(os.path.join(DATA, "02_checkm2", "mag_summary.tsv"), sep='\t')
    fam_counts = df.groupby(['family', 'cohort']).size().unstack(fill_value=0)
    fam_counts['total'] = fam_counts.sum(axis=1)
    fam_counts = fam_counts.sort_values('total', ascending=False).head(15)
    fam_counts = fam_counts.drop(columns='total')
    fig, ax = plt.subplots(figsize=(10, 7))
    fam_counts.plot(kind='barh', stacked=True, ax=ax,
                    color=[CRC_COLOR if c == 'CRC' else HEALTHY_COLOR for c in fam_counts.columns])
    ax.invert_yaxis()
    ax.set_xlabel("Number of MAGs"); ax.set_ylabel("Bacterial family")
    ax.set_title("Figure 5. Top 15 bacterial families among 226 MAGs\n(GTDB-Tk r232, 99.6% species-level resolution)",
                 fontsize=13, fontweight='bold')
    ax.legend(title='Cohort')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "fig5_gtdbtk_taxonomy.png"), dpi=SAVE_DPI, bbox_inches='tight')
    plt.savefig(os.path.join(OUT, "fig5_gtdbtk_taxonomy.pdf"), bbox_inches='tight')
    plt.close()
    print("Saved fig5_gtdbtk_taxonomy.{png,pdf}")

if __name__ == "__main__":
    print(f"Output folder: {OUT}\n")
    print("Generating Figure 1: Jaccard concordance...")
    fig1_jaccard()
    print("Generating Figure 2: Per-metabolite concordance...")
    fig2_metabolite_concordance()
    print("Generating Figure 3: Spearman heatmap...")
    fig3_spearman_heatmap()
    print("Generating Figure 4: CheckM2 quality...")
    fig4_checkm2_quality()
    print("Generating Figure 5: GTDB-Tk taxonomy...")
    fig5_gtdbtk_taxonomy()
    print(f"\nAll figures written to: {OUT}")
