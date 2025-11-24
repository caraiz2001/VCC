#%% Load filtered counts
import pandas as pd
import numpy as np
#%%
# Load the filtered counts data. it is a compressed numpy object
filtered_counts_path = 'data/filtered_counts.npz'
try:
    filtered_counts = np.load(filtered_counts_path, allow_pickle=True)
    print("Filtered counts data loaded successfully.")
    # Load the filtered counts data into a DataFrame only if loading was successful
    if all(key in filtered_counts for key in ['target', 'gene_names', 'target_sample_names']):
        target_counts_df = pd.DataFrame(
            filtered_counts['target'],
            columns=filtered_counts['gene_names'],
            index=filtered_counts['target_sample_names']
        )
        control_counts_df = pd.DataFrame(
            filtered_counts['control'],
            columns=filtered_counts['gene_names'],
            index=filtered_counts['control_sample_names']
        )
        print("Target counts DataFrame shape:", target_counts_df.shape)
        print("Control counts DataFrame shape:", control_counts_df.shape)
    else:
        print("Error: One or more required keys are missing in the filtered_counts file.")
except FileNotFoundError:
    print(f"Error: The file {filtered_counts_path} was not found.")

# %%
# I want to save control + target counts to a csv file called all_counts.csv
all_counts_df = pd.concat([control_counts_df, target_counts_df])
all_counts_df.to_csv('data/all_counts.csv')

# %%
# I want generate a samplesheet witht the metadata for each sample in all_counts_df
# The samplesheet should have the following columns:
# sample_name, condition (control or target), gene (the gene that was perturbed, or 'control' for control samples)
sample_names = all_counts_df.index.tolist()
conditions = ['control'] * len(control_counts_df) + ['target'] * len(target_counts_df)
target_genes = ['non-targeting'] * len(control_counts_df) + ['ANTXR1'] * len(target_counts_df)  # Replace 'perturbed_gene' with actual gene names if available

# %%
# I want to save the samplesheet to a csv file called samplesheet.csv
samplesheet_df = pd.DataFrame({
    'sample_name': sample_names,
    'condition': conditions,
    'gene': target_genes
})
samplesheet_df.to_csv('data/samplesheet.csv', index=False)

# %%
