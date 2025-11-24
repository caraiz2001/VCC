#%% LOAD LIBRARIES
import pyperclip as clip
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import networkx as nx

#%%
# LOAD DATA 

# File with N DE genes per targeted gene
all_genes = pd.read_csv('data/de_metrics_results.csv', sep=',')
# File with protein-protein interaction database
pp_interaction_data = pd.read_csv('data/9606.protein.links.v12.0.txt', sep=' ')
# File with Ensemble protein ID and gene name
mapping_ids = pd.read_csv('data/biomart_report.txt', sep='\t')

# %%
# Extract the gene names from all_genes (column name: target_gene)
all_gene_names = all_genes['target_gene'].tolist()

# copy all_gene_names to clipboard
clip.copy('\n'.join(all_gene_names))

# For each gene, load the corresponding DE data
de_data = {}
for gene in all_gene_names:
    try:
        de_data[gene] = pd.read_csv(f'data/de_results_per_gene/{gene}_de_genes.tsv', sep='\t')
    except FileNotFoundError:
        print(f'File for gene {gene} not found.')


# # Remove the '9606.' prefix from the 'protein1' and 'protein2' columns
# pp_interaction_data['protein1'] = pp_interaction_data['protein1'].str.replace('9606.', '', regex=False)
# pp_interaction_data['protein2'] = pp_interaction_data['protein2'].str.replace('9606.', '', regex=False)

# # Overwrite the original file with the modified data
# pp_interaction_data.to_csv('data/9606.protein.links.v12.0.txt', sep=' ', index=False)

# Remove non-unique rows from the mapping_ids DataFrame
mapping_ids = mapping_ids.drop_duplicates(subset=['Protein stable ID', 'Gene name'])
# %%
# unique number of proteins
unique_proteins = pd.concat([pp_interaction_data['protein1'], pp_interaction_data['protein2']]).unique()
num_unique_proteins = len(unique_proteins)

# %%
# Look at the distribution of scores 
plt.figure(figsize=(10, 6))
plt.hist(pp_interaction_data['combined_score'], bins=100, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of Protein-Protein Interaction Scores')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('plots/ppi_score_distribution.png', dpi=300, bbox_inches='tight')

# %%
# Dictionary to save the interactions for each gene
gene_interactions_dict = {}
compare_distributions = pd.DataFrame(columns=['gene_name', 'num_DE', 'num_DE_interactions', 'num_non_DE_interactions'])
# Iterate over all genes in all_gene_names
for gene in all_gene_names:
    # Get the corresponding protein ID for the current gene
    protein_id = mapping_ids[mapping_ids['Gene name'] == gene]['Protein stable ID'].values
    if len(protein_id) == 0:
        print(f"Protein ID not found for gene: {gene}")
        continue
    protein_id = protein_id[0]

    # for prot_id in protein_id:
    #     print(len(pp_interaction_data[pp_interaction_data['protein1'] == prot_id]))

    # Subset the data to only include interactions with the current gene's protein ID
    gene_interactions = pp_interaction_data[pp_interaction_data['protein1'] == protein_id]

    # # Add a column with the gene name of the second protein in the interaction
    # gene_interactions = gene_interactions.merge(mapping_ids[['Protein stable ID', 'Gene name']],
    #                                             left_on='protein2', right_on='Protein stable ID', how='left')
    
    # Save this to a file
    #gene_interactions.to_csv(f'data/pp_interactions/{gene}_interactions.tsv', sep='\t', index=False)

    print(f"Number of interactions for gene {gene}: {len(gene_interactions)}")
    
    gene_interactions_dict[gene] = {
        'num_interactions': len(gene_interactions),}

    # Add 2 columns with the gene names from mapping_ids dataframe (gene_name_1 and gene_name_2)
    gene_interactions = gene_interactions.merge(mapping_ids[['Protein stable ID','Gene name']],
                                                left_on='protein1', right_on='Protein stable ID', how='left')
    gene_interactions.rename(columns={'Gene name': 'gene_name_1'}, inplace=True)

    # Remove the 'Protein stable ID' column as it's no longer needed
    gene_interactions.drop(columns=['Protein stable ID'], inplace=True)

    gene_interactions = gene_interactions.merge(mapping_ids[['Protein stable ID','Gene name']],
                                                left_on='protein2', right_on='Protein stable ID', how='left')
    gene_interactions.rename(columns={'Gene name': 'gene_name_2'}, inplace=True)
    
    # Remove the 'Protein stable ID' column as it's no longer needed
    gene_interactions.drop(columns=['Protein stable ID'], inplace=True)

    # Load DE genes for the current gene
    try:
        de_genes_file = f'data/de_results_per_gene/{gene}_de_genes.tsv'
        gene_de_genes = pd.read_csv(de_genes_file, sep='\t')
    except FileNotFoundError:
        print(f"DE genes file not found for gene: {gene}")
        continue

    # Merge the DE genes with the interactions
    gene_interactions_merged = gene_de_genes.merge(gene_interactions[['gene_name_2', 'combined_score']],
                                                   left_on='feature', right_on='gene_name_2')

    # Plot statistic against combined score
    plt.figure(figsize=(10, 6))
    plt.scatter(gene_interactions_merged['statistic'], gene_interactions_merged['combined_score'],
                alpha=0.5, edgecolors='w', s=100)
    plt.title(f'{gene} DE Genes vs Protein-Protein Interaction Score')
    plt.xlabel('Statistic (DE Genes)')
    plt.ylabel('Protein-Protein Interaction Score')
    plt.grid()
    #plt.show()
    plt.savefig(f'plots/scatter/{gene}_de_genes_vs_ppi_score.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Add FDR category
    gene_interactions_merged['fdr_category'] = gene_interactions_merged['fdr'].apply(lambda x: 'fdr < 0.05' if x < 0.05 else 'fdr >= 0.05')

    # Boxplot of interaction scores
    plt.figure(figsize=(10, 6))
    gene_interactions_merged = gene_interactions_merged.dropna(subset=['combined_score'])

    fdr_less_005 = gene_interactions_merged[gene_interactions_merged['fdr_category'] == 'fdr < 0.05']['combined_score']
    n_less = len(fdr_less_005)
    fdr_greater_equal_005 = gene_interactions_merged[gene_interactions_merged['fdr_category'] == 'fdr >= 0.05']['combined_score']
    n_more = len(fdr_greater_equal_005)

    label1 = f'DE genes (n={n_less})'
    label2 = f'non-DE genes (n={n_more})'

    # From the all_genes dataframe, get the number of DE genes for the current gene
    n_de_genes = int(all_genes[all_genes['target_gene'] == gene]['number_of_genes'].values[0])

    if not fdr_less_005.empty and not fdr_greater_equal_005.empty:
        plt.boxplot([fdr_less_005, fdr_greater_equal_005],
                    labels=[label1, label2],)
        plt.title(f'Protein-Protein Interaction Scores for {gene} ({n_de_genes}) DE Genes')
        plt.ylabel('Protein-Protein Interaction Score')
        plt.grid()
        #plt.show()
        plt.savefig(f'plots/ppscores/{gene}_de_genes_ppi_score_boxplot.png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        print(f"One or both subsets are empty for gene: {gene}. Skipping boxplot.")


    # Do a statistic test to see if the interaction scores are significantly different between DE genes and non-DE genes
    
    if not fdr_less_005.empty and not fdr_greater_equal_005.empty:
        t_stat, p_value = stats.ttest_ind(fdr_less_005, fdr_greater_equal_005, equal_var=False)
        print(f'T-test for {gene}: t-statistic = {t_stat}, p-value = {p_value}')

        # Save the results to a df 
        new_row = pd.DataFrame([{
            'gene_name': gene,
            'num_DE': n_de_genes,
            'num_DE_interactions': n_less,
            'num_non_DE_interactions': n_more,
            't_statistic': t_stat,
            'p_value': p_value
        }])

        compare_distributions = pd.concat([compare_distributions, new_row], ignore_index=True)

#%%
# Sort the compare_distributions dataframe by num_DE
compare_distributions = compare_distributions.sort_values(by='num_DE', ascending=False)    

# Merge the compare_distributions df with the mapping_ids (Gene description)
compare_distributions = compare_distributions.merge(mapping_ids[['Gene name', 'Gene description']],
                                                    left_on='gene_name', right_on='Gene name',
                                                    how='left', suffixes=('', '_desc'))

# Sort by p_value
compare_distributions = compare_distributions.sort_values(by='p_value', ascending=True)

# Remove duplicates based on 'gene_name'
compare_distributions = compare_distributions.drop_duplicates(subset=['gene_name'])

# Save the compare_distributions dataframe to a file
compare_distributions.to_csv('data/compare_distributions_ppi_devsnonde.tsv', sep='\t', index=False)

# Keep only the first 8 columns
compare_distributions = compare_distributions.iloc[:, :8]
# load the manual annotations

# %%
# -------- LET'S INSPECT SYMMETRIC INTERACTIONS --------
# for each entry of 
de_genes_per_gene = {}
for gene in de_data.keys():
    if gene in de_data:
        de_genes = de_data[gene]
        de_genes_filtered = de_genes[de_genes['fdr'] < 0.05]
        de_genes_per_gene[gene] = de_genes_filtered['feature'].tolist()
        print(f'Gene: {gene}, Number of DE genes with FDR < 0.05: {len(de_genes_filtered)}')
    else:
        print(f'No DE data for gene: {gene}')


# %%
# For each gene in all_gene_names (targeted genes):
# 1. Get the DE genes for the gene from de_genes_per_gene
# 2. Filter for those that are appear in all_gene_names (de_gene_targeted)
# 3. Get the DE genes for de_gene_targeted from de_genes_per_gene
# 4. See if gene is in the DE genes of de_gene_targeted
# 5. If it is, sum 1 to the count for that gene
# 6. Save the results to a dataframe with the following columns: 'Target Gene', 'Number of DE genes that are as well targets', 'Number of symmetric relationships'

symmetric_counts = {gene: 0 for gene in all_gene_names}
total_counts = 0
sim_counts = 0
de_genes_as_target_counts = 0
de_genes_non_target_counts = 0
de_genes_total_counts = 0
non_simmetric_counts = 0
simmetric_and_interaction_counts = 0
simmetric_no_interaction_counts = 0

for gene in all_gene_names:
    if gene in de_genes_per_gene:
        de_genes = de_genes_per_gene[gene]
        # Filter for DE genes that are also in all_gene_names
        de_genes_targeted = [g for g in de_genes if g in all_gene_names]
        print(f'Gene: {gene}, Number of DE that appear as target: {len(de_genes_targeted)}')
        de_genes_as_target_counts += len(de_genes_targeted)
        de_genes_non_target_counts += (len(de_genes) - len(de_genes_targeted))
        de_genes_total_counts += len(de_genes)
        
        # For each DE gene that is a target as well, check if the original gene is in its DE genes
        for de_gene_targ in de_genes_targeted:
            if de_gene_targ in de_genes_per_gene and gene in de_genes_per_gene[de_gene_targ]:
                #symmetric_counts[gene] += 1
                sim_counts += 1
                # If the de_gene_targ has a interaction score > 600 with the original gene, count is incremented
                interaction_scores = pp_interaction_data[(pp_interaction_data['gene_name_1'] == gene) & 
                                                         (pp_interaction_data['gene_name_2'] == de_gene_targ)]['combined_score'].values
                if len(interaction_scores) > 0 and interaction_scores[0] > 600:
                    simmetric_and_interaction_counts += 1
                else:
                    simmetric_no_interaction_counts += 0
                
                total_counts += 1
            else:
                total_counts += 1
                non_simmetric_counts += 1
print(f'Total number of genes: {total_counts}')
print(f'Total number of symmetric relationships: {sim_counts}')
print(f'percentage of symmetric relationships: {sim_counts / total_counts * 100:.2f}%')
print(f'Total number of DE genes that are also targets: {de_genes_as_target_counts}')
print(f'Total number of DE genes that are not targets: {de_genes_non_target_counts}')
print(f'Total number of DE genes: {de_genes_total_counts}')
print(f'Total number of symmetric DE relationships with interaction score > 600: {simmetric_and_interaction_counts}')
print(f'Total number of symmetric DE relationships without interaction score > 600: {simmetric_no_interaction_counts}')
print(f'Total number of non-symmetric DE relationships: {non_simmetric_counts}')


# %%%
# Add 2 columns to pp_interaction data with the gene names from mapping_ids dataframe (gene_name_1 and gene_name_2)
pp_interaction_data = pp_interaction_data.merge(mapping_ids[['Protein stable ID', 'Gene name']],
                                                left_on='protein1', right_on='Protein stable ID', how='left')
pp_interaction_data.rename(columns={'Gene name': 'gene_name_1'}, inplace=True)
# Remove the 'Protein stable ID' column as it's no longer needed
pp_interaction_data.drop(columns=['Protein stable ID'], inplace=True)
pp_interaction_data = pp_interaction_data.merge(mapping_ids[['Protein stable ID', 'Gene name']],
                                                left_on='protein2', right_on='Protein stable ID', how='left')
pp_interaction_data.rename(columns={'Gene name': 'gene_name_2'}, inplace=True)
# Remove the 'Protein stable ID' column as it's no longer needed
pp_interaction_data.drop(columns=['Protein stable ID'], inplace=True)



#%% SIMILAR, BUT CONSIDERING THE INTERACTIONS
# For each gene in all_gene_names (targeted genes):
# 1. Get the protein ID from the mapping file
# 2. Get the interactions from the pp_interaction_data
# 3. Filter for significant interactions (combined_score > 600)
# 4. Save them name of the significant interactions to a dictionary with the target gene as key
interactions_per_gene = {}
total_significant_interactions = 0
for gene in all_gene_names:
    gene_interactions = pp_interaction_data[pp_interaction_data['gene_name_1'] == gene]
    significant_interactions = gene_interactions[gene_interactions['combined_score'] > 600]
    interactions_per_gene[gene] = significant_interactions['gene_name_2'].tolist()
    total_significant_interactions += len(significant_interactions)

#%%
import networkx as nx
from collections import deque

#%%
# Step 1: Build an undirected graph from PPI data
G = nx.Graph()
filtered_ppI_df = pp_interaction_data[pp_interaction_data['combined_score'] >= 600]
for _, row in filtered_ppI_df.iterrows():
    G.add_edge(row['gene_name_1'], row['gene_name_2'], weight=row['combined_score'])

# Step 2: BFS function to find nodes within max_hops from each target
def bfs_up_to_k_hops(graph, source, max_hops):
    visited = set([source])
    queue = deque([(source, 0)])
    result = set()
    
    while queue:
        node, depth = queue.popleft()
        if depth >= max_hops:
            continue
        for neighbor in graph.neighbors(node):
            if neighbor not in visited:
                visited.add(neighbor)
                result.add(neighbor)
                queue.append((neighbor, depth + 1))
    return result

valid_targets = [t for t in all_gene_names if t in G.nodes]
indirect_interactors_dict = {
    target: bfs_up_to_k_hops(G, target, max_hops=2)
    for target in valid_targets
}



# Result
print("Indirect interactors within 10 hops per CRISPR target:")
for target, interactors in indirect_interactors_dict.items():
    print(f"{target} --> {interactors}")



# %%
# For each gene in all_gene_names (targeted genes) create a dictionary with the genes that are both in the interactions and in the DE genes
interactions_and_de_genes = {}
size = 0
counts_total_de_genes = 0
counts_total_significant_interactions = 0
counts_sig_and_de_genes = 0

for gene in all_gene_names:
    if gene in de_genes_per_gene:
        # get DE genes for the gene
        de_genes = de_genes_per_gene[gene]
        # get interaction genes for the gene
        interactions = interactions_per_gene[gene]
        
        counts_total_de_genes += len(de_genes)
        counts_total_significant_interactions += len(interactions)

        # Find the intersection between the DE genes and the interactions
        common_genes = set(de_genes_targeted).intersection(set(interactions))
        counts_sig_and_de_genes += len(common_genes)

        print(f'Gene: {gene}, Number of DE genes that are also in significant interactions: {len(common_genes)}')
        interactions_and_de_genes[gene] = list(common_genes)

    else:
        print(f'No DE data for gene: {gene}')

# %%
nuclear_genes = [
    "HMGN1", "SALL4", "SIX4", "SAFB", "HIRA", "PRDM14", "ZNF714", "ZNF426", "ZNF562",
    "EID2", "SSBP1", "DHX36", "HMGB2", "METTL17", "XRCC4", "TADA1", "EWSR1", "POLB",
    "METTL3", "MED12", "IRF3", "PHF14", "KLF10", "FOXH1", "MTA1", "IKBKG", "HAT1", "RRM1",
    "TET1", "CENPO", "SUPT4H1", "DNMT1", "TCF3", "KDM1A", "ETV4", "SMARCA4", "MED1",
    "MED13", "TAF13", "STAT6", "USP22", "MED24", "PHF10", "STAT1", "ATM", "SOX2", "KDM2B",
    "PBX1", "ZNF593", "RNF2", "PMS1", "ARID1A"
]
de_and_interact_counts = 0
not_de_and_interact_counts = 0
total_inspected = 0
si_and_de_butnontargeted = 0
# inspect simmetric relationship
for gene in interactions_and_de_genes.keys():
    # if len(de_genes_per_gene[gene]) >= 1000:
    #     print(f'Gene: {gene} has too many DE genes ({len(de_genes_per_gene[gene])}) to inspect symmetric relationships.')
    #     continue
    # if gene not in nuclear_genes:
    #     continue
    print(f'Inspecting gene: {gene} with {len(interactions_and_de_genes[gene])} interacting DE genes out of {len(de_genes_per_gene[gene])} DE genes.')
    if len(interactions_and_de_genes[gene]) > 0:
        for related_gene in interactions_and_de_genes[gene]:
            total_inspected += 1
            print(f'Inspecting symmetric relationship between {gene} and {related_gene}')
            if related_gene in all_gene_names:
                if related_gene in interactions_and_de_genes and gene in interactions_and_de_genes[related_gene]:
                    print(f'Gene: {gene} has a symmetric relationship with {related_gene}.')
                    de_and_interact_counts += 1
                else:
                    print(f'Gene: {gene} has no symmetric relationship with {related_gene}.')
                    not_de_and_interact_counts += 1
            else:
                print(f'Gene: {related_gene} is not a targeted gene.')
                si_and_de_butnontargeted += 1
    else:
        print(f'Gene: {gene} has no interactions.')
        #not_de_and_interact_counts += 1
        # Uncomment the following line if you want to print genes with no symmetric relationships
        #print(f'Gene: {gene} has no symmetric relationships.')
    #else:
        #print(f'Gene: {gene} has no symmetric relationships.')
print(f'Total number of symmetric relationships: {de_and_interact_counts}')
print(f'Total number of non-symmetric relationships: {not_de_and_interact_counts}')
print(f'percentage of symmetric relationships: {de_and_interact_counts / (de_and_interact_counts + not_de_and_interact_counts) * 100:.2f}%')

# %%
# filter mapped_ids to only include the genes that are in all_gene_names
mapped_ids_filtered = mapping_ids[mapping_ids['Gene name'].isin(all_gene_names)]
# remove duplicates based on 'Gene name'
mapped_ids_filtered = mapped_ids_filtered.drop_duplicates(subset=['Gene name'])
# Save the filtered mapping_ids to a file
mapped_ids_filtered.to_csv('data/mapped_ids_filtered.tsv', sep='\t', index=False)
# %%

# Similar analysis but a bit more specific
# I want to create a dataframe with the following columns:
# - Target Gene: gene name
# - Other Gene: gene name of the other gene in the pair
# - other_DE: boolean indicating if the other gene is a DE gene of the target gene
# - interaction_score: interaction score between the target gene and the other gene
# - target_DE: boolean indicating if the target gene is a DE gene of the other gene
# - symmetric: boolean indicating if the relationship is symmetric

# In order to do this:
# 1. Create a df with the columns mentioned above
# 2. Write all gene names in the 'Target Gene' column

#results_df = pd.DataFrame(columns=['Target Gene', 'Other Gene', 'num_DE', 'num_DE_also_target', 'is_O_DE_of_T', 'interaction_score', 'is_T_DE_of_O', 'symmetric'])

results_df = pd.DataFrame(columns=['Target Gene', 'Other Gene', 'num_DE', 'num_DE_also_target',
                                  'is_O_DE_of_T', 'is_T_DE_of_O', 'symmetric'])
# Write all pairs (150 target genes * 150 other genes) to the df
for target_gene in all_gene_names:
    print(f'Processing target gene: {target_gene}')
    # Get the DE genes for the target gene
    if target_gene in de_genes_per_gene:
        de_genes = de_genes_per_gene[target_gene]
    else:
        de_genes = []
        print(f'No DE data for target gene: {target_gene}')
    
    # # Get the interactions for the target gene
    # if target_gene in interactions_per_gene:
    #     interactions = interactions_per_gene[target_gene]
    # else:
    #     interactions = []
    
    # Iterate over all other genes
    for other_gene in all_gene_names:
        if other_gene == target_gene:
            continue  # Skip self-comparison
        
        # Check if the other gene is a DE gene of the target gene
        is_O_DE_of_T = other_gene in de_genes
        
        # Check if the target gene is a DE gene of the other gene
        is_T_DE_of_O = target_gene in de_genes_per_gene.get(other_gene, [])
        
        # Check if the relationship is symmetric
        symmetric = is_O_DE_of_T and is_T_DE_of_O

        # Count DE genes that are also targeted
        num_DE_also_targeted = len([g for g in de_genes if g in all_gene_names])
        
        new_row = pd.DataFrame([{
            'Target Gene': target_gene,
            'Other Gene': other_gene,
            'num_DE': len(de_genes),
            'num_DE_also_target': num_DE_also_targeted,
            'is_O_DE_of_T': is_O_DE_of_T,
            #'interaction_score': interaction_score,
            'is_T_DE_of_O': is_T_DE_of_O,
            'symmetric': symmetric
        }])
        results_df = pd.concat([results_df, new_row], ignore_index=True)


#%%
# Remove the rows gene_name_1, gene_name_2, and combined_score from the results_df
results_df = results_df.drop(columns=['gene_name_1', 'gene_name_2', 'combined_score'], errors='ignore')

#%%
# Add a column with the protein IDs for the target gene and the other gene
results_df_mapped = results_df.merge(mapping_ids[['Gene name', 'Protein stable ID']],
                              left_on='Target Gene', right_on='Gene name',
                              how='left', suffixes=('', '_target'))

results_df_mapped = results_df_mapped.merge(mapping_ids[['Gene name', 'Protein stable ID']],
                              left_on='Other Gene', right_on='Gene name',
                              how='left', suffixes=('', '_other'))

# Merge the results_df with the pp_interaction_data to get the interaction score
results_df_mapped = results_df_mapped.merge(pp_interaction_data[['protein1','protein2', 'combined_score']],
                              left_on=['Protein stable ID', 'Protein stable ID_other'],
                              right_on=['protein1', 'protein2'],
                              how='left')


# Add a c
# %%
# Check how many interaction scores are NaN
nan_interaction_scores = results_df_mapped['combined_score_y'].isna().sum()

# Print the first 10 rows that are not NaN
print("First 10 rows with non-NaN interaction scores:")
print(results_df[results_df['combined_score'].notna()].head(10))
# %%

results_noNaN = results_df_mapped[results_df_mapped['combined_score'].notna()]

# %%
# print the rows in which target gene is CHMP3 and other gene is AKT2
print(results_noNaN[(results_noNaN['Target Gene'] == 'CHMP3') & (results_noNaN['Other Gene'] == 'AKT2')])


# If other is in de of target (results_noNaN['is_O_DE_of_T'] == True) and combined_score > 600, then print the row
print("Rows where Other Gene is DE of Target Gene and interaction score > 600:")
print(results_noNaN[(results_noNaN['is_O_DE_of_T'] == True) & (results_noNaN['combined_score'] > 600)])

#%%
# Define significance thresholds to test
sig_thresholds = list(range(0, int(results_noNaN['combined_score'].max()), 50))
  # Set a significance threshold for the combined score

# create empty df to store results
empty_df = pd.DataFrame(columns=['Significance Threshold', 'Number of Symmetric Rows', 'Percentage of Symmetric Rows'])
# Try with different significance threshold (from 0 to max, in steps of 50)
for t in sig_thresholds:
    symmetric_rows = results_noNaN[(results_noNaN['is_O_DE_of_T'] == True) & 
                               (results_noNaN['is_T_DE_of_O'] == True) &
                               (results_noNaN['combined_score'] > t)]
    print(f"Number of symmetric rows: {len(symmetric_rows)}")
    percentage_symmetric = len(symmetric_rows) / len(results_noNaN[(results_noNaN['is_O_DE_of_T'] == True) & (results_noNaN['combined_score'] > t)])* 100
    print(f"Percentage of symmetric rows: {percentage_symmetric:.2f}%")
    # Append the results to the empty_df
    new_row = pd.DataFrame([{
        'Significance Threshold': t,
        'Number of Symmetric Rows': len(symmetric_rows),
        'Percentage of Symmetric Rows': percentage_symmetric
    }])
    empty_df = pd.concat([empty_df, new_row], ignore_index=True)


#%%

# Plot the significance threshold vs percentage of symmetric rows with dual y-axes
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot percentage of symmetric rows on the first y-axis
ax1.plot(empty_df['Significance Threshold'], empty_df['Percentage of Symmetric Rows'], marker='o', linestyle='-', color='blue', label='Percentage of Symmetric Relationships')
ax1.set_xlabel('Significance Threshold')
ax1.set_ylabel('Percentage of Symmetric Relationships (%)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.grid()

# Create a second y-axis for the number of symmetric rows
ax2 = ax1.twinx()
ax2.plot(empty_df['Significance Threshold'], empty_df['Number of Symmetric Rows'], marker='o', linestyle='-', color='red', label='Number of Symmetric Relationships')
ax2.set_ylabel('Number of Symmetric Relationships', color='red')
ax2.tick_params(axis='y', labelcolor='red')

# Add a title and legend
fig.suptitle('Percentage and Number of Symmetric Relationships vs Significance Threshold')
fig.tight_layout()
plt.savefig('plots/symmetric_rows_vs_significance_threshold.png', dpi=300, bbox_inches='tight')


# %%
import pandas as pd


# %%
# load the manual annotations
manual_annotations = pd.read_csv('data/All_Annotated_Proteins___150_.csv', sep=',')
# %%
# load compare distributions
compare_distributions_df = pd.read_csv('data/compare_distributions_ppi_devsnonde.tsv', sep='\t')
# %%



# ---------------- PROPAGATION TESTING ----------------
# The idea is that we want to test if the propagation of DE genes through the PPI network is significant.
# That is. If Gene B is a DE gene of Gene A, and Gene B is physicaly interacting with Gene C, then we want to see if Gene C is also a DE gene.

# For this, we will:
# 1. For each gene in all_gene_names, get the genes that are indirectly interacting with it (i.e. genes that are up to 10 hops away).
# 2. For each of those genes, check if they are DE genes of the original gene (= precision)
# 3. percentage of DE genes that are indirectly interacting with the original gene (= recall)
#%% 
# Create a dictionary to store the propagation results
propagation_results = {}
# Iterate over all genes in all_gene_names
for gene in all_gene_names:
    print(f'Processing gene: {gene}')
    # Get the interactions for the gene from interactions_per_gene
    if gene in indirect_interactors_dict:
        interactions = indirect_interactors_dict[gene]
    else:
        print(f'No interactions for gene: {gene}')
        continue

    # Get the DE genes for the gene from de_genes_per_gene
    if gene in de_genes_per_gene:
        de_genes = de_genes_per_gene[gene]
    else:
        print(f'No DE genes for gene: {gene}')
        de_genes = []
        continue

    # Calculate precision and recall
    precision = len(set(de_genes).intersection(interactions)) / len(interactions) if interactions else 0
    recall = len(set(de_genes).intersection(interactions)) / len(de_genes) if de_genes else 0
    
    # Save the results to the dictionary
    propagation_results[gene] = {
        'precision': precision,
        'recall': recall,
        'num_indirect_interactions': len(interactions),
        'num_de_genes': len(de_genes)
    }
# Convert the results to a DataFrame
propagation_results_df = pd.DataFrame.from_dict(propagation_results, orient='index')


# %%
# For each target gene, take the DE genes
# For each DE gene, see how many of the direct interactions are also DE genes
# Save the results to a dictionary with the target gene as key

are_de_close = {}
for gene in all_gene_names:
    print(f'Processing gene: {gene}')
    # Get the DE genes for the gene from de_genes_per_gene
    if gene in de_genes_per_gene:
        de_genes = de_genes_per_gene[gene]
    else:
        print(f'No DE genes for gene: {gene}')
        de_genes = []
        continue
    
    num_de_close = 0
    total_interactions = 0
    for de_gene in de_genes:
        print(f'Processing DE gene: {de_gene}')
        # Get the interactions for the DE gene from interactions_per_gene
        if de_gene in interactions_per_gene:
            interactions = interactions_per_gene[de_gene]
        else:
            print(f'No interactions for DE gene: {de_gene}')
            continue
        # Calculate how many of the interactions are DE genes
        num_de_close += len(set(de_genes).intersection(interactions))
        total_interactions += len(interactions)
    
    # Save the results to the dictionary
    are_de_close[gene] = {
        'num_de_close': num_de_close,
        'num_interactions': len(interactions),
        'num_de_genes': len(de_genes),
        'total_interactions': total_interactions,
    }
# %%
# convert the results to a DataFrame
are_de_close_df = pd.DataFrame.from_dict(are_de_close, orient='index')
# Compute the percentage of DE genes that are close
are_de_close_df['percentage_de_close'] = (are_de_close_df['num_de_close'] / 
                                           are_de_close_df['total_interactions']) * 100
print(are_de_close_df.head())
# %%
# Histogram of the percentage of DE genes that are close
plt.figure(figsize=(10, 6))
plt.hist(are_de_close_df['percentage_de_close'], bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of Percentage of DE Genes That Are Close')
plt.xlabel('Percentage of DE Genes That Are Close (%)')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('plots/percentage_de_genes_close_histogram.png', dpi=300, bbox_inches='tight')

# %%
# Do a df with the following columns:
# - Target Gene
# - Number of DE genes
# - Number of DE genes that are close
# - Percentage of DE genes that are close
# - Number of interactions (taken from the interactions_per_gene dictionary)

propagation_summary_df = pd.DataFrame(columns=['Target Gene', 'Number of DE Genes',
                                              'Number of DE Genes That Are Close',
                                              'Percentage of DE Genes That Are Close',
                                              'Number of Interactions'])
for gene, data in propagation_results.items():
    #print(data['num_de_genes'])
    new_row = pd.DataFrame([{
        'Target Gene': gene,
        'Number of DE Genes': data['num_de_genes'],
        'Number of DE Genes That Are Close': data['num_indirect_interactions'],
        'Percentage of DE Genes That Are Close': data['recall'] * 100,  # recall is already a percentage
        'Number of Interactions': data['num_indirect_interactions']
    }])
    propagation_summary_df = pd.concat([propagation_summary_df, new_row], ignore_index=True)

propagation_summary_df = propagation_summary_df.rank()

propagation_summary_df.corr('spearman')

# %%
# Scatter plot of the number of DE genes vs the number of interactions
plt.figure(figsize=(10, 6))
plt.scatter(propagation_summary_df['Number of DE Genes'], propagation_summary_df['Number of Interactions'],
            alpha=0.5, edgecolors='w', s=100)
plt.title('Number of DE Genes vs Number of Interactions')
plt.xlabel('Number of DE Genes')
plt.ylabel('Number of Interactions')
plt.grid()
plt.savefig('plots/de_genes_vs_interactions_scatter.png', dpi=300, bbox_inches='tight')


# %%
# Rank the df by number of DE genes
propagation_summary_df = propagation_summary_df.sort_values(by='Number of DE Genes', ascending=False)

# Select the bottom 30% of the rows
bottom_30_percent = propagation_summary_df.head(int(len(propagation_summary_df) * 0.3))

# Compute correlation between the number of DE genes and the number of interactions
correlation = bottom_30_percent['Number of DE Genes'].corr(bottom_30_percent['Number of Interactions'], method='spearman')
print(f'Correlation between the number of DE genes and the number of interactions: {correlation:.2f}')

# %%
