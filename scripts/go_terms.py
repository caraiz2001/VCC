#%%
import mygene
import pandas as pd
import numpy as np

#%%
mg = mygene.MyGeneInfo()

# Load data
crispr_data = pd.read_csv('data/crispr_effect_results.csv')
de_data = pd.read_csv('data/de_metrics_results.csv')

# Merge data on gene names
merged_data = pd.merge(de_data, crispr_data, on='target_gene', how='inner')


# List of your gene symbols
genes = merged_data['target_gene'].tolist()

#%%
# Query MyGene.info
results = mg.querymany(genes, scopes='symbol', fields='go', species='human')

# Initialize dictionaries
go_bp = {}
go_mf = {}
go_cc = {}

for res in results:
    if 'notfound' in res and res['notfound']:
        continue

    gene = res['query']
    go = res.get('go', {})

    # Helper to extract terms
    def extract_terms(entry):
        if isinstance(entry, list):
            return [e['term'] for e in entry[:3]]
        else:
            return [entry['term']]

    bp_terms = extract_terms(go.get('BP', [])) if 'BP' in go else []
    mf_terms = extract_terms(go.get('MF', [])) if 'MF' in go else []
    cc_terms = extract_terms(go.get('CC', [])) if 'CC' in go else []

    go_bp[gene] = "; ".join(bp_terms)
    go_mf[gene] = "; ".join(mf_terms)
    go_cc[gene] = "; ".join(cc_terms)


# %%
merged_data['GO_BP'] = merged_data['target_gene'].map(go_bp)
merged_data['GO_MF'] = merged_data['target_gene'].map(go_mf)
merged_data['GO_CC'] = merged_data['target_gene'].map(go_cc)


# %%
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

# Your list of target genes
genes = merged_data['target_gene'].tolist()

# Query g:Profiler
# 'GO:BP', 'GO:MF', 'GO:CC' will be mapped to GO slim by default
results = gp.profile(
    organism='hsapiens',
    query=genes,
    sources=['GO:BP', 'GO:MF', 'GO:CC'],
    no_evidences=True,
    user_threshold=1.0,
    all_results=True,
    domain_scope='annotated',
    slim=True
)

#%%
# Create dictionaries for each GO aspect
go_slim_bp = (
    results[results['source'] == 'GO:BP']
    .groupby('query')['name']
    .apply(lambda x: "; ".join(x))
    .to_dict()
)

go_slim_mf = (
    results[results['source'] == 'GO:MF']
    .groupby('query')['name']
    .apply(lambda x: "; ".join(x))
    .to_dict()
)

go_slim_cc = (
    results[results['source'] == 'GO:CC']
    .groupby('query')['name']
    .apply(lambda x: "; ".join(x))
    .to_dict()
)


# %%
# Optional: Filter to only GO slim terms
# g:Profiler auto-maps to a slimmed set, but you can filter further if needed

# Group all slim terms per gene
slim_map = results.groupby('query')['name'].apply(lambda terms: "; ".join(terms)).to_dict()

# Add to your DataFrame
merged_data['GO_slim'] = merged_data['target_gene'].map(slim_map)

# %%
