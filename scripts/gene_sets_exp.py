"""
We would like to identify which Gene Pairs are more likely to be DE together.
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%%
# Load the DE data for all genes
de_file_names = pd.read_csv('data/de_results_per_gene/de_results_filename.txt', header=None)
# Load all_gene_names.txt
all_gene_names = pd.read_csv('data/all_genes_names.txt', header=None)[0].tolist()

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

de_gene_lists = {}
for gene, df in de_data.items():
    print(f"Processing gene: {gene}")
    significant_de = df[df['fdr'] < 0.05]
    total_de = len(significant_de)
    print(f"Gene: {gene}, Total DE genes with FDR < 0.05: {total_de}")
    if total_de > 0:
        de_gene_lists[gene] = significant_de['feature'].tolist()
    else:
        de_gene_lists[gene] = []

# %%
# I want to create a df where the rows are all genes (from all_gene_names.txt) and the columns are targeted genes (from de_gene_lists keys)
# The values are 1 if the gene is DE when the targeted gene is perturbed, 0 otherwise
all_genes_df = pd.DataFrame(0, index=all_gene_names, columns=de_gene_lists.keys())
for target_gene, de_genes in de_gene_lists.items():
    for gene in de_genes:
        if gene in all_genes_df.index:
            all_genes_df.at[gene, target_gene] = 1
# Save the df to a csv file
#all_genes_df.to_csv('data/gene_pair_de_matrix.csv')

# Add a column with the sum of each row
all_genes_df['total_de_count'] = all_genes_df.sum(axis=1)
# Sort by the total_de_count column
all_genes_df_sorted = all_genes_df.sort_values(by='total_de_count', ascending=False)
# Save the sorted df to a csv file

# Histogram of the total_de_count column
plt.figure(figsize=(10, 6))
plt.hist(all_genes_df_sorted['total_de_count'], bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Distribution of DE Counts Across All Genes')
plt.xlabel('Number of Times Gene is DE')
plt.ylabel('Number of Genes')
plt.grid()
plt.show()

all_genes_matrix = all_genes_df_sorted.drop(columns=['total_de_count'])

# THis is jaccard index between experiments
from sklearn.metrics import jaccard_score

n = all_genes_matrix.shape[1]
jaccard_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i, n):
        score = jaccard_score(all_genes_matrix.iloc[:, i], all_genes_matrix.iloc[:, j])
        jaccard_matrix[i, j] = jaccard_matrix[j, i] = score

# Now for jaccard index between genes
import pandas as pd
from collections import defaultdict
from itertools import combinations

df = all_genes_matrix

# df: genes x experiments, binary
gene_sets = {gene: set(df.columns[df.loc[gene] == 1]) for gene in df.index}


# Inverted index: experiment -> genes DE
experiment_to_genes = defaultdict(set)
for gene, exps in gene_sets.items():
    for e in exps:
        experiment_to_genes[e].add(gene)

# Candidate pairs: only genes that co-occur at least once
pairs = defaultdict(int)
for genes in experiment_to_genes.values():
    for g1, g2 in combinations(genes, 2):
        pairs[(g1, g2)] += 1

THRESHOLD = 0.1  # keep only stronger co-occurrences
jaccard_scores = {}
for (g1, g2), inter in pairs.items():
    union = len(gene_sets[g1]) + len(gene_sets[g2]) - inter
    score = inter / union
    if score >= THRESHOLD:
        jaccard_scores[(g1, g2)] = score

jaccard_df = pd.DataFrame(
    [(g1, g2, score) for (g1, g2), score in jaccard_scores.items()],
    columns=["Gene1", "Gene2", "Jaccard"]
)

# median jaccard score
median_jaccard = jaccard_df['Jaccard'].median()
print(f"Median Jaccard Score: {median_jaccard}")
mean_jaccard = jaccard_df['Jaccard'].mean()
print(f"Mean Jaccard Score: {mean_jaccard}")

# FIlter for jaccard > 0.25
jaccard_df_filtered = jaccard_df[jaccard_df['Jaccard'] > 0.25]

# Sort by jaccard score
jaccard_df_sorted = jaccard_df_filtered.sort_values(by='Jaccard', ascending=False)



# --------------------------
# Weired code efficient testing
# --------------------------

import os
import math
import gc
import uuid
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def compute_fisher_pvalue_from_counts(overlap, count_a, count_b, total_experiments):
    """One-sided Fisher exact test for enrichment of overlap."""
    a = int(overlap)
    b = int(count_a - overlap)
    c = int(count_b - overlap)
    d = int(total_experiments - (count_a + count_b - overlap))
    # Protect against degenerate cells
    if a < 0 or b < 0 or c < 0 or d < 0:
        return 1.0
    _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(p_value)

def stream_significant_pairs_pass1(
    binary_gene_by_experiment,                       # DataFrame: genes x experiments, values in {0,1}
    output_dir,
    block_size=512,
    min_experiments_per_gene=5,
    min_shared_experiments=2,
    optional_keep_top_k_per_gene=None,               # e.g., 100; if None, keep all that pass p filter
    optional_max_pairs_per_chunk=2_000_000,          # for big outputs, split chunks
    optional_max_pvalue_to_keep=0.05                 # early prune before FDR to cut disk size
):
    """
    Pass 1:
    - Filters genes by activity
    - Computes overlaps in blocks without building the full square matrix
    - Writes candidate pairs with p-values to one or more parquet/csv chunk files
    Returns a list of file paths of the chunks written.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df = binary_gene_by_experiment.astype(np.int8)
    total_experiments = df.shape[1]

    # Filter by gene frequency
    per_gene_counts = df.sum(axis=1)
    valid_genes = per_gene_counts[per_gene_counts >= min_experiments_per_gene].index
    df = df.loc[valid_genes]
    per_gene_counts = per_gene_counts.loc[valid_genes]

    # Precompute a transposed int8 matrix for fast dot products
    # shape: experiments x genes
    matrix_experiments_by_genes = df.T.values  # int8
    gene_names = df.index.to_numpy()
    gene_count_array = per_gene_counts.to_numpy()

    n_genes = df.shape[0]
    chunk_paths = []
    rows_buffer = []
    rows_in_buffer = 0

    # Iterate in blocks over the first axis (rows)
    for start in range(0, n_genes, block_size):
        stop = min(start + block_size, n_genes)
        block_idx = np.arange(start, stop)

        # block: genes_in_block x experiments
        block = df.iloc[block_idx].values  # int8
        # overlaps: (genes_in_block x experiments) @ (experiments x genes) => genes_in_block x n_genes
        # result entries are integer overlaps 0..total_experiments
        overlaps = block @ matrix_experiments_by_genes  # int32

        # For each gene i in the block, consider only j>i to avoid duplicates
        for local_i, i in enumerate(block_idx):
            # vector of overlaps to all genes
            ov = overlaps[local_i]
            # enforce j > i
            j_indices = np.arange(n_genes)
            mask_upper = j_indices > i
            # apply minimum overlap and early p-value pruning later
            candidate_js = j_indices[mask_upper & (ov >= min_shared_experiments)]
            if candidate_js.size == 0:
                continue

            overlaps_i = ov[candidate_js]
            count_i = gene_count_array[i]
            counts_j = gene_count_array[candidate_js]

            # compute p values for all candidate pairs for this i
            # vectorized loop for clarity; you can numba if needed
            pvals = np.fromiter(
                (compute_fisher_pvalue_from_counts(int(overlap_ij), int(count_i), int(count_j), int(total_experiments))
                 for overlap_ij, count_j in zip(overlaps_i, counts_j)),
                dtype=float,
                count=candidate_js.size
            )

            # early prune to reduce disk writes
            keep_mask = pvals <= optional_max_pvalue_to_keep
            if optional_keep_top_k_per_gene is not None and keep_mask.any():
                # keep top-k smallest p values for this gene
                idx_sorted = np.argsort(pvals[keep_mask])
                idx_sorted = idx_sorted[:optional_keep_top_k_per_gene]
                keep_indices = np.flatnonzero(keep_mask)[idx_sorted]
            else:
                keep_indices = np.flatnonzero(keep_mask)

            if keep_indices.size == 0:
                continue

            for k in keep_indices:
                j = candidate_js[k]
                rows_buffer.append((
                    gene_names[i], gene_names[j],
                    int(overlaps_i[k]),
                    int(count_i),
                    int(counts_j[k]),
                    float(pvals[k])
                ))
                rows_in_buffer += 1

            # spill to disk if buffer large
            if rows_in_buffer >= optional_max_pairs_per_chunk:
                chunk_df = pd.DataFrame(rows_buffer, columns=[
                    "gene_one", "gene_two",
                    "shared_experiments",
                    "count_gene_one", "count_gene_two",
                    "p_value"
                ])
                chunk_path = Path(output_dir) / f"pairs_pass1_{uuid.uuid4().hex}.parquet"
                try:
                    chunk_df.to_parquet(chunk_path, index=False)
                except Exception:
                    # fallback if parquet not available
                    chunk_path = Path(output_dir) / f"pairs_pass1_{uuid.uuid4().hex}.csv.gz"
                    chunk_df.to_csv(chunk_path, index=False, compression="gzip")
                chunk_paths.append(str(chunk_path))
                rows_buffer.clear()
                rows_in_buffer = 0

        # keep memory tight
        del overlaps, block
        gc.collect()

    # flush remainder
    if rows_buffer:
        chunk_df = pd.DataFrame(rows_buffer, columns=[
            "gene_one", "gene_two",
            "shared_experiments",
            "count_gene_one", "count_gene_two",
            "p_value"
        ])
        chunk_path = Path(output_dir) / f"pairs_pass1_{uuid.uuid4().hex}.parquet"
        try:
            chunk_df.to_parquet(chunk_path, index=False)
        except Exception:
            chunk_path = Path(output_dir) / f"pairs_pass1_{uuid.uuid4().hex}.csv.gz"
            chunk_df.to_csv(chunk_path, index=False, compression="gzip")
        chunk_paths.append(str(chunk_path))
        rows_buffer.clear()

    return chunk_paths

def apply_global_false_discovery_rate_pass2(
    chunk_paths,
    output_final_path,
    false_discovery_rate=0.05,
    optional_min_shared_experiments=None
):
    """
    Pass 2:
    - Loads all chunks’ p values
    - Applies one global Benjamini–Hochberg correction
    - Writes final significant pairs
    """
    # Load only needed columns first
    frames = []
    for pth in chunk_paths:
        if pth.endswith(".parquet"):
            dfc = pd.read_parquet(pth, columns=["p_value"])
        else:
            dfc = pd.read_csv(pth, usecols=["p_value"])
        frames.append(dfc)
    all_p = pd.concat(frames, ignore_index=True)["p_value"].to_numpy()
    if all_p.size == 0:
        pd.DataFrame(columns=[
            "gene_one", "gene_two",
            "shared_experiments",
            "count_gene_one", "count_gene_two",
            "p_value", "q_value"
        ]).to_parquet(output_final_path, index=False)
        return output_final_path

    # Global BH
    _, qvals, _, _ = multipletests(all_p, alpha=false_discovery_rate, method="fdr_bh")
    # We need to reattach q values to the original rows without loading everything at once
    writer = None
    q_offset = 0
    out_ext = Path(output_final_path).suffix.lower()

    for pth in chunk_paths:
        # read full chunk
        if pth.endswith(".parquet"):
            dfc = pd.read_parquet(pth)
        else:
            dfc = pd.read_csv(pth)
        n = len(dfc)
        dfc["q_value"] = qvals[q_offset:q_offset + n]
        q_offset += n

        # optional extra filter
        if optional_min_shared_experiments is not None:
            dfc = dfc[dfc["shared_experiments"] >= optional_min_shared_experiments]

        dfc = dfc[dfc["q_value"] <= false_discovery_rate]

        if dfc.empty:
            continue

        # append to final file
        if out_ext == ".parquet":
            if writer is None:
                dfc.to_parquet(output_final_path, index=False)
                writer = "parquet"
            else:
                # append by reading existing and concatenating chunk at a time is costly;
                # better to collect then write once if you expect many chunks to survive.
                # Simpler approach: write csv.gz instead for true append.
                tmp = pd.read_parquet(output_final_path)
                pd.concat([tmp, dfc], ignore_index=True).to_parquet(output_final_path, index=False)
        else:
            mode = "w" if writer is None else "a"
            header = writer is None
            dfc.to_csv(output_final_path, index=False, mode=mode, header=header)
            writer = "csv"

    return output_final_path





# %%
# df is your genes x experiments binary dataframe
chunk_paths = stream_significant_pairs_pass1(
    df,
    output_dir="cooccurrence_chunks",
    block_size=512,
    min_experiments_per_gene=5,
    min_shared_experiments=2,
    optional_keep_top_k_per_gene=100,         # keeps the 100 smallest p values per gene, change as needed
    optional_max_pairs_per_chunk=1_000_000,   # adjust to your disk and patience
    optional_max_pvalue_to_keep=0.05
)

final_path = apply_global_false_discovery_rate_pass2(
    chunk_paths,
    output_final_path="significant_gene_pairs.parquet",
    false_discovery_rate=0.05,
    optional_min_shared_experiments=2
)
print("Final file:", final_path)

# %%
