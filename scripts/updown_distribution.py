# Load the file with the DE file names
import pandas as pd

#%%
de_file_names = pd.read_csv('data/de_results_per_gene/de_results_filename.txt', header=None)
# %%
# Each row contains a file name, so we can use it to load the DE data
de_data = {}
for file_name in de_file_names[0]:
    try:
        gene_name = file_name.split('_')[0]  # Extract the gene name from the file name
        de_data[gene_name] = pd.read_csv(f'data/de_results_per_gene/{file_name}', sep='\t')
        print(f"Loaded DE data for {gene_name}.")
    except FileNotFoundError:
        print(f"File {file_name} not found.")
    except Exception as e:
        print(f"An error occurred while loading {file_name}: {e}")
# %%
# For each gene:
#  1.count the number of DE genes with fdr < 0.05
#  2. from the significant, count the number of upregulated and downregulated genes (log2fc > 0 and log2fc < 0)
#  3. Compute the percentage of upregulated/totalDE and downregulated/totalDE genes
de_summary = {}
for gene, df in de_data.items():
    print(f"Processing gene: {gene}")
    significant_de = df[df['fdr'] < 0.05]
    total_de = len(significant_de)
    print(f"Gene: {gene}, Total DE genes with FDR < 0.05: {total_de}")
    if total_de > 0:
        upregulated = significant_de[significant_de['fold_change'] > 1]
        downregulated = significant_de[significant_de['fold_change'] <= 1]

        upreg_count = len(upregulated)
        downreg_count = len(downregulated)
        
        upreg_percentage = (upreg_count / total_de) * 100
        downreg_percentage = (downreg_count / total_de) * 100
        
        de_summary[gene] = {
            'total_DE': total_de,
            'upregulated': upreg_count,
            'downregulated': downreg_count,
            'upreg_percentage': upreg_percentage,
            'downreg_percentage': downreg_percentage
        }
        #break
    else:
        de_summary[gene] = {
            'total_DE': 0,
            'upregulated': 0,
            'downregulated': 0,
            'upreg_percentage': 0.0,
            'downreg_percentage': 0.0
        }
        #break
#%%
# Convert the summary to a DataFrame for easier analysis
de_summary_df = pd.DataFrame.from_dict(de_summary, orient='index')
de_summary_df.reset_index(inplace=True)
de_summary_df.rename(columns={'index': 'gene'}, inplace=True)

# Sort the Dataframe by percentage of upregulated genes in descending oreder
de_summary_df.sort_values(by='upreg_percentage', ascending=False, inplace=True)
# Save the summary to a CSV file
de_summary_df.to_csv('data/de_distribution_summary.tsv', sep='\t', index=False)

# %%

# Print the row of df in which feature is == gene
def print_gene_row(df, gene):
    row = df[df['feature'] == gene]
    if not row.empty:
        print(row)
    else:
        print(f"No data found for gene: {gene}")
print_gene_row(df, gene)

# %%
# Plot a histogram of upreg_percentage
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.hist(de_summary_df['upreg_percentage'], bins=30, color='blue', alpha=0.7)
plt.title('Distribution of Upregulated Genes Percentage')
plt.xlabel('Upregulated Percentage')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.show()

# %%
