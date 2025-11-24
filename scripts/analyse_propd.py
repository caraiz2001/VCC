import pandas as pd

# Load MED12_counts.csv from the data directory
#%%
data_path = 'data/MED12_counts.csv'
try:
    med12_counts = pd.read_csv(data_path)
    print("MED12 counts data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {data_path} was not found.")
# %%
print(med12_counts.head())
# print shape of the dataframe
print(f"Data shape: {med12_counts.shape}")
# %%
# Print minimum and maximum of each column

# Select only numeric columns for min/max calculations
numeric_df = med12_counts.select_dtypes(include='number')
min_values = numeric_df.min()
max_values = numeric_df.max()

# Print summary of the distribution of the min values across columns (average of the min values)
print(f"Average of minimum values across columns: {min_values.mean()}")
# Print minimum of the minimum values across columns
print(f"Minimum of minimum values across columns: {min_values.min()}")
# Print summary of the distribution of the max values across columns (average of the max values)
print(f"Average of maximum values across columns: {max_values.mean()}")
# Print minimum of the maximum values across columns
print(f"Minimum of maximum values across columns: {max_values.min()}")

# %%
# Load the gsea_propd.csv file
gsea_path = 'data/gsea_propd_nobug.csv'
all_gene_names_path = 'data/all_genes_names.txt'
ANTXR1_de_path = 'data/ANTXR1_de.tsv'
weighted_connectivity_path = 'data/weighted_con_ANTXR1_nobug.csv'
try:
    gsea_propd = pd.read_csv(gsea_path)
    print("GSEA propd data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {gsea_path} was not found.")

# Load the all_gene_names.csv file
try:
    all_gene_names = pd.read_csv(all_gene_names_path)
    print("All gene names data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {all_gene_names_path} was not found.")

try:
    antrx1_de = pd.read_csv(ANTXR1_de_path, sep='\t')
    print("ANTXR1 DE data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {ANTXR1_de_path} was not found.")

# Load the weighted connectivity matrix

try:
    weighted_connectivity = pd.read_csv(weighted_connectivity_path)
    print("Weighted connectivity data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {weighted_connectivity_path} was not found.")

#%%
# Add a column named propd_index to all_gene_names that starts from 1 and increments by 1 for each row
all_gene_names['propd_index'] = range(1, len(all_gene_names) + 1)
# Sort gsea_propd by the gene column (index)
gsea_propd_sorted = gsea_propd.sort_values(by='gene')

# Merge with all_gene_names. The index of all_gene_names matches the 'gene' column in gsea_propd_sorted
merged_data = pd.merge(gsea_propd_sorted, all_gene_names, left_on='gene', right_on='propd_index', how='left')

# Merge with the weighted connectivity matrix
merged_data = pd.merge(merged_data, weighted_connectivity, left_on='gene', right_on='gene_index', how='left')

# Filter the antrx1_de data to include only significantly DE genes (fdr < 0.05)
significant_antrx1_de = antrx1_de[antrx1_de['fdr'] < 0.05]

# Print the number of significant DE genes
print(f"Number of significant DE genes in ANTXR1: {significant_antrx1_de.shape[0]}")


# Add a column to the merged_data indicating if the gene is significantly DE in ANTXR1
merged_data['is_sig_de'] = merged_data['gene_name'].isin(significant_antrx1_de['feature'])

# Sort data by the is_sig_de column to have True values first
#merged_data_sortedtrue = merged_data.sort_values(by=['is_sig_de', 'ES'], ascending=False)

# Sort by ES_pos
merged_data_sorted = merged_data.sort_values(by='ES_pos', ascending=False)


# I want to track the number of True Positives (is_sig_de == True) when sorted by ES_pos. Plot the cumulative sum of True Positives against the index of the sorted DataFrame.
import matplotlib.pyplot as plt
# Calculate cumulative sum of True Positives
cumulative_true_positives = merged_data_sorted['is_sig_de'].cumsum()
# Reset the index to use the index of the sorted DataFrame for plotting
merged_data_sorted.reset_index(drop=True, inplace=True)
# Plot the cumulative sum of True Positives
plt.figure(figsize=(10, 6))
plt.plot(merged_data_sorted.index, cumulative_true_positives, linestyle='-', color='blue')
plt.title('Cumulative Sum of True Positives (is_sig_de == True) Sorted by ES_pos')
plt.xlabel('Index of Sorted DataFrame')
plt.ylabel('Cumulative Sum of True Positives')
plt.grid()
plt.show()


# Same plot but sorting for weighted connectivity
merged_data_sorted_weighted = merged_data.sort_values(by='weighted_connectivity', ascending=False)
# Calculate cumulative sum of True Positives
cumulative_true_positives_weighted = merged_data_sorted_weighted['is_sig_de'].cumsum()
# Reset the index to use the index of the sorted DataFrame for plotting
merged_data_sorted_weighted.reset_index(drop=True, inplace=True)
# Plot the cumulative sum of True Positives
plt.figure(figsize=(10, 6))
plt.plot(merged_data_sorted_weighted.index, cumulative_true_positives_weighted, linestyle='-', color='green')
plt.plot(merged_data_sorted.index, cumulative_true_positives, linestyle='-', color='purple')
plt.title('Agreement of DE genes between PROPD and Wilcoxson ')
plt.xlabel('Genes sorted by propd DE-genewise score')
plt.ylabel('Cumulative Sum of True DE Genes (FDR < 0.05 according to Wilcoxson)')
# Add legend of colors
plt.legend(['Weighted Connectivity', 'ES_pos'])
# X axis should go from 0 to 2000
#plt.xlim(0, 2000)
plt.grid()
plt.show()

# Read compressed file
file = 'data/gsea_propd_nobug.csv.gz'


