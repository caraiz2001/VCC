# Load AKT2 data
import pandas as pd
import numpy as np

#%%
# Load data
akt2_de = pd.read_csv('data/akt2_de.tsv', sep = '\t')
# %%
# Count the number of genes with fdr < 0.05
num_de_genes = akt2_de[akt2_de['fdr'] < 0.05].shape[0]
print(f'Number of DE genes with FDR < 0.05: {num_de_genes}')

# %% Save the column 'feature' to a file
akt2_de['feature'].to_csv('data/all_genes_names.txt', index=False, header=False)

# %%
# print the name of the the top 10 genes with the lowest fdr
top_de_genes = akt2_de.sort_values(by='fdr').head(10)
print("Top 10 DE genes with lowest FDR:")
print(top_de_genes[['feature']])
# %%
# load the file with all genes 
all_genes = pd.read_csv('data/de_metrics_results.csv', sep=',')
# %%
# Extract the gene names from all_genes (column name: target_gene)
all_gene_names = all_genes['target_gene'].tolist()

# For each gene, load the corresponding DE data
de_data = {}
for gene in all_gene_names:
    try:
        de_data[gene] = pd.read_csv(f'data/de_results_per_gene/{gene}_de_genes.tsv', sep='\t')
    except FileNotFoundError:
        print(f'File for gene {gene} not found.')
#%%
# Sort all genes by the number of DE genes
all_genes_sorted = all_genes.sort_values(by='number_of_genes', ascending=False)
# Didive in 3 groups: top 33%, middle 33%, bottom 33%
top_33 = all_genes_sorted.head(int(len(all_genes_sorted) * 0.33))
middle_33 = all_genes_sorted.iloc[int(len(all_genes_sorted) * 0.33):int(len(all_genes_sorted) * 0.66)]
bottom_33 = all_genes_sorted.tail(int(len(all_genes_sorted) * 0.33))
# %%
# For the bottom 33% of genes:
# 1. Extrac the de data from the de_data dictionary
# 2. Count the number of DE genes for each gene (fdr < 0.05)
# 3. Print the number of DE genes for each gene
# 4. Save them to a dictionary with the gene name as key and the the list of DE genes as value

bottom_33_de_data = {}
for gene in bottom_33['target_gene']:
    if gene in de_data:
        de_genes = de_data[gene]
        de_genes_filtered = de_genes[de_genes['fdr'] < 0.05]
        bottom_33_de_data[gene] = de_genes_filtered['feature'].tolist()
        print(f'Gene: {gene}, Number of DE genes with FDR < 0.05: {len(de_genes_filtered)}')
    else:
        print(f'No DE data for gene: {gene}')
# %%
# Calculate how many times each gene appears differentially expressed in the bottom 33% of genes
from collections import Counter
gene_counts = Counter()
for de_genes in bottom_33_de_data.values():
    gene_counts.update(de_genes)  

# Do a df with gene names and their counts
gene_counts_df = pd.DataFrame(gene_counts.items(), columns=['gene', 'count'])
# Sort by count
gene_counts_df = gene_counts_df.sort_values(by='count', ascending=False)
# %%
# inspect the distribution of counts
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.hist(gene_counts_df['count'], bins=30, alpha=0.7, edgecolor='black')
plt.title('Distribution of DE Gene Counts in Bottom 33%')
plt.xlabel('How many times the gene appears as DE in the 50 genes with the lowest number of DE genes')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('plots/de_gene_counts_distribution.png', dpi=300, bbox_inches='tight')

# %%
# Same for the top 33% of genes
top_33_de_data = {}
for gene in top_33['target_gene']:
    if gene in de_data:
        de_genes = de_data[gene]
        de_genes_filtered = de_genes[de_genes['fdr'] < 0.05]
        top_33_de_data[gene] = de_genes_filtered['feature'].tolist()
        print(f'Gene: {gene}, Number of DE genes with FDR < 0.05: {len(de_genes_filtered)}')
    else:
        print(f'No DE data for gene: {gene}')


# %%
# Calculate how many times each gene appears differentially expressed in the top 33% of genes
top_gene_counts = Counter()
for de_genes in top_33_de_data.values():
    top_gene_counts.update(de_genes)

# %%
# Do a df with gene names and their counts and plot
top_gene_counts_df = pd.DataFrame(top_gene_counts.items(), columns=['gene', 'count'])
# Sort by count
top_gene_counts_df = top_gene_counts_df.sort_values(by='count', ascending=False)    

# Inspect the distribution of counts
plt.figure(figsize=(10, 6))
plt.hist(top_gene_counts_df['count'], bins=30, alpha=0.7, edgecolor='black')
plt.title('Distribution of DE Gene Counts in Top 33%')
plt.xlabel('How many times the gene appears as DE in the 50 genes with the highest number of DE genes')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('plots/top_de_gene_counts_distribution.png', dpi=300, bbox_inches='tight')
# %%
# Same for the middle 33% of genes
middle_33_de_data = {}
for gene in middle_33['target_gene']:
    if gene in de_data:
        de_genes = de_data[gene]
        de_genes_filtered = de_genes[de_genes['fdr'] < 0.05]
        middle_33_de_data[gene] = de_genes_filtered['feature'].tolist()
        print(f'Gene: {gene}, Number of DE genes with FDR < 0.05: {len(de_genes_filtered)}')
    else:
        print(f'No DE data for gene: {gene}')

# %%
# Calculate how many times each gene appears differentially expressed in the middle 33% of genes
middle_gene_counts = Counter()
for de_genes in middle_33_de_data.values():
    middle_gene_counts.update(de_genes)
# Do a df with gene names and their counts
middle_gene_counts_df = pd.DataFrame(middle_gene_counts.items(), columns=['gene', 'count'])
# Sort by count
middle_gene_counts_df = middle_gene_counts_df.sort_values(by='count', ascending=False)
# Inspect the distribution of counts
plt.figure(figsize=(10, 6))
plt.hist(middle_gene_counts_df['count'], bins=30, alpha=0.7, edgecolor='black')
plt.title('Distribution of DE Gene Counts in Middle 33%')
plt.xlabel('How many times the gene appears as DE in the 50 genes with the middle')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('plots/middle_de_gene_counts_distribution.png', dpi=300, bbox_inches='tight')
# %%
# for the gene AKT2, print the name of the top 10 DE genes with the lowest fdr
interesting_gene = 'METTL14'
if interesting_gene in de_data:
    akt2_de = de_data[interesting_gene]
akt2_de_sorted = akt2_de.sort_values(by='fdr').head(10)
print("Top 10 DE genes for AKT2 with lowest FDR:")
cb = akt2_de_sorted[['feature']]
import pyperclip as clip
clip.copy(cb.to_string(index=False))
      
# %%
