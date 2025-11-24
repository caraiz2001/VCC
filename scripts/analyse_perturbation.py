# Load perturbation analysis data of MED12 gene
import pandas as pd

#%%
# Load the perturbation data
perturbation_data_path = 'data/MED12_Xpert_all_genes.csv'
try:
    perturbation_data = pd.read_csv(perturbation_data_path)
    print("Perturbation analysis data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {perturbation_data_path} was not found.")
# %%
# Transpose the data to have genes as rows and samples as columns. Rows are genes, columns are samples
perturbation_data_t = perturbation_data.set_index('Unnamed: 0').T

# Change the 'Unnamed: 0' column to 'gene_name'
perturbation_data_t.index.name = 'gene_name'
perturbation_data_t.reset_index(inplace=True)

# %%
# See if there are some 'conserved' genes across samples.
# Built a different df indicating the number of samples in which each gene has a positive value
summary_df = perturbation_data_t.set_index('gene_name').apply(lambda x: x > 0).sum(axis=1).reset_index()
summary_df.columns = ['gene_name', 'num_samples_positive']

# Add a column with the same but negative values
summary_df['num_samples_negative'] = perturbation_data_t.set_index('gene_name').apply(lambda x: x < 0).sum(axis=1).values

# Add a column with the number of samples with zero values
summary_df['num_samples_zero'] = perturbation_data_t.set_index('gene_name').apply(lambda x: x == 0).sum(axis=1).values

# Add a column with the percentage of samples with positive values
summary_df['percentage_positive'] = (summary_df['num_samples_positive'] / perturbation_data_t.shape[1]) * 100
# Add a column with the percentage of samples with negative values
summary_df['percentage_negative'] = (summary_df['num_samples_negative'] / perturbation_data_t.shape[1]) * 100

# Add a column with the difference between positive and negative percentages (absolute value)
summary_df['percentage_difference_abs'] = (summary_df['percentage_positive'] - summary_df['percentage_negative']).abs()

summary_df['percentage_difference'] = (summary_df['percentage_positive'] - summary_df['percentage_negative'])

# Sort the summary_df by the absolute value of the percentage difference abs
summary_df_sorted = summary_df.sort_values(by='percentage_difference_abs', ascending=False)

# Remove those genes that have zero values in all samples
summary_df_sorted_filtered = summary_df_sorted[summary_df_sorted['num_samples_zero'] != 2766]
# %%
# Plot histogram of the percentage difference (absolute value)
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.hist(summary_df_sorted['percentage_difference_abs'], bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of Percentage Difference (Absolute Value)')
plt.xlabel('Percentage Difference (Absolute Value)')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.axvline(x=0, color='red', linestyle='--', label='Zero Difference')
plt.legend()
plt.show()
# %%
# Compute the covariance matrix
# Load the MED12 counts
med12_counts_path = 'data/MED12_counts.csv'
try:
    med12_counts = pd.read_csv(med12_counts_path)
    print("MED12 counts data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {med12_counts_path} was not found.")

# Set the first column as index
med12_counts.set_index('Unnamed: 0', inplace=True)

#%%
# Load covariance data
cov_data_path = 'data/MED12_cov_vector.csv'
try:
    cov_data = pd.read_csv(cov_data_path)
    print("Covariance data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {cov_data_path} was not found.")

#%%
# Load all genes names
all_gene_names_path = 'data/all_genes_names.txt'
try:
    all_gene_names = pd.read_csv(all_gene_names_path, sep='\t')
    print("All gene names loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {all_gene_names_path} was not found.")
# %%
# Merge the covariance data with all gene names (set gene_name as index)
all_gene_names.set_index('gene_name', inplace=True)
# Set index of all_gene_names as index of cov_data
cov_data.index = all_gene_names.index
cov_data.index.name = 'gene_name'

# Print the covariance with the gene_name MED12
print("Covariance for MED12:")
print(cov_data.loc['MED12'])

# %%
# Merge the summary_df with the cov_data according to the gene_name
summary_df_sorted_filtered.set_index('gene_name', inplace=True)
merged_summary_cov = summary_df_sorted_filtered.join(cov_data, how='inner', rsuffix='_cov')


# Add a column indicating if the sign in percentage_difference matches the sign in the covariance
merged_summary_cov['sign_match'] = (merged_summary_cov['percentage_difference'] > 0) == (merged_summary_cov['MED12'] > 0)


# Sort the merged_summary_cov by the absolute value of the percentage difference
merged_summary_cov_sorted = merged_summary_cov.sort_values(by='percentage_difference_abs', ascending=False)

# How many genes have a sign match?
num_sign_match = merged_summary_cov_sorted['sign_match'].sum()
print(f"Number of genes with sign match: {num_sign_match}")

# Save the merged summary with covariance to a CSV file
merged_summary_cov_sorted.to_csv('data/merged_summary_cov_sorted.csv')

# Add a column with the cumulative sum of the sign matches
merged_summary_cov_sorted['cumulative_sign_match'] = merged_summary_cov_sorted['sign_match'].cumsum()
# Add a column with the index (from 1 to n) of the sorted DataFrame
merged_summary_cov_sorted['index_num'] = range(1, len(merged_summary_cov_sorted) + 1)
#%% Plot the cumulative sum of sign matches against the index of the sorted DataFrame
plt.figure(figsize=(10, 6))
plt.plot(merged_summary_cov_sorted['index_num'], merged_summary_cov_sorted['cumulative_sign_match'], marker='o', linestyle='-', color='blue')
plt.title('Cumulative Sum of Sign Matches (Percentage Difference Sign Matches Covariance Sign)')
plt.xlabel('Gene Index (Sorted by Percentage Difference Absolute Value)')
plt.ylabel('Cumulative Sum of Sign Matches')
plt.grid()
plt.show()
# %%

# Define 'conserved' genes as those with percentage_difference_abs > 60% and non conserved as those with percentage_difference_abs < 20%
merged_summary_cov_sorted['conserved'] = merged_summary_cov_sorted['percentage_difference_abs'] > 60
merged_summary_cov_sorted['non_conserved'] = merged_summary_cov_sorted['percentage_difference_abs'] < 20

# Count the number of conserved and non-conserved genes
num_conserved = merged_summary_cov_sorted['conserved'].sum()
num_non_conserved = merged_summary_cov_sorted['non_conserved'].sum()
print(f"Number of conserved genes: {num_conserved}")
print(f"Number of non-conserved genes: {num_non_conserved}")

# Compute the percentage of 'sign_match' == True for conserved and non-conserved genes separetly
conserved_sign_match_percentage = merged_summary_cov_sorted[merged_summary_cov_sorted['conserved']]['sign_match'].mean() * 100
non_conserved_sign_match_percentage = merged_summary_cov_sorted[merged_summary_cov_sorted['non_conserved']]['sign_match'].mean() * 100

print(f"Percentage of sign matches in conserved genes: {conserved_sign_match_percentage:.2f}%")
print(f"Percentage of sign matches in non-conserved genes: {non_conserved_sign_match_percentage:.2f}%")



#%% Load MED12 differntial expression data
med12_de_path = 'data/de_results_per_gene/MED12_de_genes.tsv'
try:
    med12_de = pd.read_csv(med12_de_path, sep='\t')
    print("MED12 DE data loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file {med12_de_path} was not found.")
# %%
# Add a column indicating if the gene is significantly DE (fdr < 0.05)
med12_de['significantly_DE'] = med12_de['fdr'] < 0.05

# Count the number of significantly DE genes
num_significantly_de_genes = med12_de['significantly_DE'].sum()
print(f"Number of significantly DE genes in MED12: {num_significantly_de_genes}")


# Merge the merged_summary_cov_sorted with the MED12 DE data (left on gene_name, right on feature).  Keep gene_name as index
med12_de.set_index('feature', inplace=True)
merged_summary_cov_de = merged_summary_cov_sorted.join(med12_de[['fdr', 'significantly_DE']], how='left', rsuffix='_de')


# print med12_de with feature  = 'MED12'
print("MED12 DE data with feature as index:")
print(med12_de[med12_de['feature'] == 'MED12'])

# Add a column with the cumulative sum of significantly DE genes
merged_summary_cov_de['cumulative_significantly_DE'] = merged_summary_cov_de['significantly_DE'].cumsum()

# Plot the cumulative sum of significantly DE genes against the index of the sorted DataFrame
plt.figure(figsize=(10, 6))
plt.plot(merged_summary_cov_de['index_num'], merged_summary_cov_de['cumulative_significantly_DE'], marker='o', linestyle='-', color='green')
plt.title('Cumulative Sum of Significantly DE Genes (fdr < 0.05)')
plt.xlabel('Gene Index (Sorted by Percentage Difference Absolute Value)')
plt.ylabel('Cumulative Sum of Significantly DE  Genes')
plt.grid()
plt.show()

# Percentage of significantly DE genes in conserved and non-conserved genes
conserved_significantly_de_percentage = merged_summary_cov_de[merged_summary_cov_de['conserved']]['significantly_DE'].mean() * 100
non_conserved_significantly_de_percentage = merged_summary_cov_de[merged_summary_cov_de['non_conserved']]['significantly_DE'].mean() * 100

print(f"Percentage of significantly DE genes in conserved genes: {conserved_significantly_de_percentage:.2f}%")
print(f"Percentage of significantly DE genes in non-conserved genes: {non_conserved_significantly_de_percentage:.2f}%")

# Percentage of sig DE overall
overall_significantly_de_percentage = merged_summary_cov_de['significantly_DE'].mean() * 100
print(f"Overall percentage of significantly DE genes: {overall_significantly_de_percentage:.2f}%")



# Do a boxplot of percentage_difference_abs for significantly DE and non-significantly DE genes
plt.figure(figsize=(10, 6))
merged_summary_cov_de.boxplot(column='percentage_difference_abs', by='significantly_DE',
                              grid=False, patch_artist=True,
                              boxprops=dict(facecolor='lightblue', color='black'),
                              medianprops=dict(color='red'),
                              whiskerprops=dict(color='black'),
                              capprops=dict(color='black'),
                              flierprops=dict(marker='o', color='black', markersize=5))
plt.title('Boxplot of Percentage Difference Absolute Value by Significantly DE Status')
plt.suptitle('')
plt.xlabel('Significantly DE Status (fdr < 0.05)')
plt.ylabel('Percentage Difference Absolute Value')
plt.xticks([1, 2], ['Not Significantly DE', 'Significantly DE'])
plt.grid(axis='y', alpha=0.75)
plt.show()


# Plot percentage_difference against covariance (column 'MED12')
plt.figure(figsize=(10, 6))
plt.scatter(merged_summary_cov_de['MED12'], merged_summary_cov_de['percentage_difference'], alpha=0.5)
plt.title('Percentage Difference vs. MED12 Covariance')
plt.xlabel('MED12 Covariance')
plt.ylabel('Percentage Difference')
plt.grid()
plt.show()

# Plot distribution of covariance values
plt.figure(figsize=(10, 6))
plt.hist(merged_summary_cov_de['MED12'], bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of MED12 Covariance Values')
plt.xlabel('MED12 Covariance')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.axvline(x=0, color='red', linestyle='--', label='Zero Covariance')
plt.legend()
plt.show()


# Print min, max, mean, and std of the covariance values
min_cov = merged_summary_cov_de['MED12_cov'].min()
max_cov = merged_summary_cov_de['MED12_cov'].max()
mean_cov = merged_summary_cov_de['MED12_cov'].mean()
std_cov = merged_summary_cov_de['MED12_cov'].std()
print(f"MED12 Covariance - Min: {min_cov}, Max: {max_cov}, Mean: {mean_cov}, Std: {std_cov}")

# Print number of positive, negative and zero covariance values
num_positive_cov = (merged_summary_cov_de['MED12_cov'] > 0).sum()
num_negative_cov = (merged_summary_cov_de['MED12_cov'] < 0).sum()
num_zero_cov = (merged_summary_cov_de['MED12_cov'] == 0).sum()
print(f"Number of positive covariance values: {num_positive_cov}")
print(f"Number of negative covariance values: {num_negative_cov}")
print(f"Number of zero covariance values: {num_zero_cov}")




# Add a column to the merged_summary_cov_de incating the rank of the covariance values (biggest value should be 1)
merged_summary_cov_de['cov_rank'] = merged_summary_cov_de['MED12'].abs().rank(ascending=False, method='min')

# Rename the columns for clarity
# percentage_difference_abs to conservation_score
merged_summary_cov_de.rename(columns={'percentage_difference_abs': 'conservation_score'}, inplace=True)
# MED12 to MED12_cov
merged_summary_cov_de.rename(columns={'MED12': 'MED12_cov'}, inplace=True)

# Add a column with the absolute value of the covariance (MED12_cov)
merged_summary_cov_de['MED12_cov_abs'] = merged_summary_cov_de['MED12_cov'].abs()


# Plot the distribution of percentage_difference (histogram )
plt.figure(figsize=(10, 6))
plt.hist(merged_summary_cov_de['percentage_difference'], bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of Percentage Difference')
plt.xlabel('Percentage Difference')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.axvline(x=0, color='red', linestyle='--', label='Zero Difference')
plt.legend()
plt.show()


# when the covariance is negative, calculate the number of negative percentage_difference values
negative_percentage_difference = merged_summary_cov_de[merged_summary_cov_de['MED12_cov'] < 0]['percentage_difference']