# Script to analyse if CRISPR and the number of DE genes are correlated

#%%
# Load libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
# set the working directory to /Users/crsitina/Documents/VCC
os.chdir('/Users/crsitina/Documents/VCC/')

# Load data
crispr_data = pd.read_csv('data/crispr_effect_results.csv')
de_data = pd.read_csv('data/de_metrics_results.csv')




# Merge data on gene names
merged_data = pd.merge(de_data, crispr_data, on='target_gene', how='inner')

# %%
# Sort by absolute value of log2 fold change
merged_data['abs_log2fc'] = merged_data['log2FC'].abs()
merged_data['log_mean_expr_control'] = np.log(merged_data['mean_expr_control'] + 1e-6)
sorted_data = merged_data.sort_values(by='abs_log2fc', ascending=False)

# %% Boxplot of the number of DE genes for top 10 and bottom 10 genes
top_10_genes = sorted_data.head(10)
bottom_10_genes = sorted_data.tail(10)

plt.figure(figsize=(10, 6))
plt.boxplot([top_10_genes['number_of_genes'], bottom_10_genes['number_of_genes']],
            labels=['Top 10 Genes with highest LFC', 'Bottom 10 Genes with lowest LFC'],)
plt.title('Number of DE Genes for Top 10 and Bottom 10 Target Genes (sorted by CRISPR effectivity)')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_boxplot.png', dpi=300, bbox_inches='tight')
# %%

# scatter plot of the number of DE genes vs abs log2 fold change
plt.figure(figsize=(10, 6))
plt.scatter(merged_data['abs_log2fc'], merged_data['number_of_genes'],
            alpha=0.5, edgecolors='w', s=100)
plt.title('Number of DE Genes vs Absolute Log2 Fold Change')
plt.xlabel('Absolute Log2 Fold Change')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_scatter.png', dpi=300, bbox_inches='tight')

#%%
# Scatter plot of number of DE genes vs abs abs_diff
plt.figure(figsize=(10, 6))
plt.scatter(merged_data['abs_diff'].abs(), merged_data['number_of_genes'],
            alpha=0.5, edgecolors='w', s=100)
plt.title('Number of DE Genes vs Absolute Difference')
plt.xlabel('Absolute Difference')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_scatter_abs_diff.png', dpi=300, bbox_inches='tight')
# %%
# scatter plot of number of DE genes vs average expression control
plt.figure(figsize=(10, 6))
plt.scatter(np.log(merged_data['mean_expr_control']), merged_data['number_of_genes'],
            alpha=0.5, edgecolors='w', s=100)
plt.title('Number of DE Genes vs Average Expression Control')
plt.xlabel('Log(Average Expression Control)')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_scatter_mean_expr_control.png', dpi=300, bbox_inches='tight')
# %%
# scatter plot of number of DE genes vs average expression crispr
plt.figure(figsize=(10, 6))
plt.scatter(np.log(merged_data['mean_expr_crispr']), merged_data['number_of_genes'],
            alpha=0.5, edgecolors='w', s=100)
plt.title('Number of DE Genes vs Average Expression CRISPR')
plt.xlabel('Log(Average Expression CRISPR)')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_scatter_mean_expr_crispr.png', dpi=300, bbox_inches='tight')
# %%
# Gene with the highest abs_diff
highest_abs_diff_gene = merged_data.loc[merged_data['abs_diff'].abs().idxmax()]
print(f"Gene with the highest absolute difference: {highest_abs_diff_gene['target_gene']} with abs_diff {highest_abs_diff_gene['abs_diff']}")
print(f"Number of DE genes: {highest_abs_diff_gene['number_of_genes']}")
print(f"Average expression control: {highest_abs_diff_gene['mean_expr_control']}")
print(f"Average expression CRISPR: {highest_abs_diff_gene['mean_expr_crispr']}")
print(f"Log2 Fold Change: {highest_abs_diff_gene['log2FC']}")
# %%
# Boxplot of number of DE genes (all genes)
plt.figure(figsize=(10, 6))
plt.boxplot(merged_data['number_of_genes'], vert=True, labels=['All Target Genes'])
plt.title('Number of DE Genes for All Target Genes')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_boxplot_all_genes.png', dpi=300, bbox_inches='tight')
# %%
# gene with the lowest number of DE genes
lowest_de_gene = merged_data.loc[merged_data['number_of_genes'].idxmin()]
print(f"Gene with the lowest number of DE genes: {lowest_de_gene['target_gene']} with {lowest_de_gene['number_of_genes']} DE genes")
# %%

# sort by number of DE genes 
sorted_by_de = merged_data.sort_values(by='number_of_genes', ascending=False)

# %%
crispr_yes = sorted_data[sorted_data['median_expr_crispr'] == 0]
crispr_no = sorted_data[sorted_data['median_expr_crispr'] != 0]

# Boxplot of number of DE genes for CRISPR yes and no
plt.figure(figsize=(10, 6))
plt.boxplot([crispr_yes['number_of_genes'], crispr_no['number_of_genes']],
            labels=['CRISPR Yes', 'CRISPR No'])
plt.title('Number of DE Genes for CRISPR Yes and No')
plt.ylabel('Number of DE Genes')
plt.grid()
plt.savefig('plots/crispr_de_boxplot_crispr_yes_no.png', dpi=300, bbox_inches='tight')
# %%
import plotly.express as px


# Assuming `merged_data` is your pandas DataFrame and includes:
# - 'mean_expr_control': average expression in control
# - 'number_of_genes': number of DE genes
# - 'target_gene': the gene name

fig = px.scatter(
    merged_data,
    x='log_mean_expr_control',
    y='number_of_genes',
    hover_name='target_gene',
    hover_data={
        'mean_expr_control': ':.2f',
        'number_of_genes': True,
        'GO_CC': True,
        'GO_MF': True,
        'GO_BP': True,
        'GO_slim': True,
        'log_mean_expr_control': False  # hide redundant value in hover
    },
    title='Number of DE Genes vs Average Expression Control (Interactive)',
    labels={
        'log_mean_expr_control': 'Log(Average Expression Control)',
        'number_of_genes': 'Number of DE Genes'
    },
    opacity=0.6,
    template='plotly_white'
)

fig.update_traces(marker=dict(size=10, line=dict(width=1, color='white')))
fig.update_layout(height=600, width=1000)
#fig.show()

# %%
fig.write_html("plots/interactive_de_vs_expr.html")
print("Plot saved to plots/interactive_de_vs_expr.html — open in browser to view it.")

# %%

