import pandas as pd
import numpy as np

# Load the datasets
p2vp0 = pd.read_csv('/Users/yitong/Documents/GitHub/ClaudeWorkshop/sc-RNAseq_data/DE.SC1.P2vP0.csv', index_col=0)
p5vp2 = pd.read_csv('/Users/yitong/Documents/GitHub/ClaudeWorkshop/sc-RNAseq_data/DE.SC1.P5vP2.csv', index_col=0)

print("="*80)
print("DIFFERENTIAL EXPRESSION OVERLAP ANALYSIS")
print("="*80)
print(f"\nDataset 1 (d0 vs d2): {len(p2vp0)} genes")
print(f"Dataset 2 (d2 vs d5): {len(p5vp2)} genes")

# Get gene lists
genes_d0vd2 = set(p2vp0.index)
genes_d2vd5 = set(p5vp2.index)

# Find overlaps
overlap = genes_d0vd2 & genes_d2vd5
print(f"\nOverlapping genes: {len(overlap)}")

# Use significance threshold
sig_threshold = 0.05

# Get significant genes with direction
sig_d0vd2_up = set(p2vp0[(p2vp0['p_val_adj'] < sig_threshold) & (p2vp0['avg_log2FC'] > 0)].index)
sig_d0vd2_down = set(p2vp0[(p2vp0['p_val_adj'] < sig_threshold) & (p2vp0['avg_log2FC'] < 0)].index)

sig_d2vd5_up = set(p5vp2[(p5vp2['p_val_adj'] < sig_threshold) & (p5vp2['avg_log2FC'] > 0)].index)
sig_d2vd5_down = set(p5vp2[(p5vp2['p_val_adj'] < sig_threshold) & (p5vp2['avg_log2FC'] < 0)].index)

print(f"\nSignificant genes (adj p < 0.05):")
print(f"  d0→d2 UP: {len(sig_d0vd2_up)} | DOWN: {len(sig_d0vd2_down)}")
print(f"  d2→d5 UP: {len(sig_d2vd5_up)} | DOWN: {len(sig_d2vd5_down)}")

# Find genes with consistent trends
consistently_up = sig_d0vd2_up & sig_d2vd5_up
consistently_down = sig_d0vd2_down & sig_d2vd5_down
reversal_up_to_down = sig_d0vd2_up & sig_d2vd5_down
reversal_down_to_up = sig_d0vd2_down & sig_d2vd5_up

print(f"\n" + "="*80)
print("CONSISTENT TREND GENES (Most Interesting)")
print("="*80)
print(f"\nConsistently UPREGULATED (d0→d2→d5): {len(consistently_up)} genes")
if consistently_up:
    result_up = []
    for gene in sorted(consistently_up, key=lambda x: p2vp0.loc[x, 'avg_log2FC'], reverse=True):
        fc_d0vd2 = p2vp0.loc[gene, 'avg_log2FC']
        fc_d2vd5 = p5vp2.loc[gene, 'avg_log2FC']
        result_up.append((gene, fc_d0vd2, fc_d2vd5))

    for gene, fc1, fc2 in result_up[:20]:  # Top 20
        print(f"  {gene:20s} | d0→d2: {fc1:+.2f} | d2→d5: {fc2:+.2f} | Total: {fc1+fc2:+.2f}")

print(f"\nConsistently DOWNREGULATED (d0→d2→d5): {len(consistently_down)} genes")
if consistently_down:
    result_down = []
    for gene in sorted(consistently_down, key=lambda x: p2vp0.loc[x, 'avg_log2FC']):
        fc_d0vd2 = p2vp0.loc[gene, 'avg_log2FC']
        fc_d2vd5 = p5vp2.loc[gene, 'avg_log2FC']
        result_down.append((gene, fc_d0vd2, fc_d2vd5))

    for gene, fc1, fc2 in result_down[:20]:  # Top 20
        print(f"  {gene:20s} | d0→d2: {fc1:+.2f} | d2→d5: {fc2:+.2f} | Total: {fc1+fc2:+.2f}")

print(f"\n" + "="*80)
print("REVERSAL TREND GENES (Direction Changes)")
print("="*80)
print(f"\nUP in d0→d2, then DOWN in d2→d5: {len(reversal_up_to_down)} genes")
if reversal_up_to_down:
    for gene in sorted(reversal_up_to_down, key=lambda x: p2vp0.loc[x, 'avg_log2FC'], reverse=True)[:10]:
        fc_d0vd2 = p2vp0.loc[gene, 'avg_log2FC']
        fc_d2vd5 = p5vp2.loc[gene, 'avg_log2FC']
        print(f"  {gene:20s} | d0→d2: {fc_d0vd2:+.2f} | d2→d5: {fc_d2vd5:+.2f}")

print(f"\nDOWN in d0→d2, then UP in d2→d5: {len(reversal_down_to_up)} genes")
if reversal_down_to_up:
    for gene in sorted(reversal_down_to_up, key=lambda x: p2vp0.loc[x, 'avg_log2FC'])[:10]:
        fc_d0vd2 = p2vp0.loc[gene, 'avg_log2FC']
        fc_d2vd5 = p5vp2.loc[gene, 'avg_log2FC']
        print(f"  {gene:20s} | d0→d2: {fc_d0vd2:+.2f} | d2→d5: {fc_d2vd5:+.2f}")

# Export results
print(f"\n" + "="*80)
print("RECOMMENDATIONS")
print("="*80)
print("\n✓ CONSISTENTLY UPREGULATED are your strongest candidates:")
print("  - These show sustained increase across the entire time course")
print("  - Likely involved in progressive differentiation/activation")
print("  - Top candidates for functional validation")

print("\n✓ CONSISTENTLY DOWNREGULATED are suppressed genes:")
print("  - May be losing function as process progresses")
print("  - Could be markers of a departing cell state")

print("\n⚠ REVERSAL genes suggest complex regulation:")
print("  - Peak expression at d2 then decline (UP→DOWN)")
print("  - Or induction at d2 after initial loss (DOWN→UP)")
print("  - May have specific temporal roles")

# Save to files for further analysis
upregulated_df = pd.DataFrame([
    {'Gene': gene, 'd0_to_d2_log2FC': p2vp0.loc[gene, 'avg_log2FC'],
     'd2_to_d5_log2FC': p5vp2.loc[gene, 'avg_log2FC'],
     'Total_log2FC': p2vp0.loc[gene, 'avg_log2FC'] + p5vp2.loc[gene, 'avg_log2FC']}
    for gene in sorted(consistently_up, key=lambda x: p2vp0.loc[x, 'avg_log2FC'], reverse=True)
])

downregulated_df = pd.DataFrame([
    {'Gene': gene, 'd0_to_d2_log2FC': p2vp0.loc[gene, 'avg_log2FC'],
     'd2_to_d5_log2FC': p5vp2.loc[gene, 'avg_log2FC'],
     'Total_log2FC': p2vp0.loc[gene, 'avg_log2FC'] + p5vp2.loc[gene, 'avg_log2FC']}
    for gene in sorted(consistently_down, key=lambda x: p2vp0.loc[x, 'avg_log2FC'])
])

upregulated_df.to_csv('/Users/yitong/Documents/GitHub/ClaudeWorkshop/sc-RNAseq_data/consistently_upregulated.csv', index=False)
downregulated_df.to_csv('/Users/yitong/Documents/GitHub/ClaudeWorkshop/sc-RNAseq_data/consistently_downregulated.csv', index=False)

print(f"\n✓ Saved results:")
print(f"  - consistently_upregulated.csv ({len(upregulated_df)} genes)")
print(f"  - consistently_downregulated.csv ({len(downregulated_df)} genes)")
