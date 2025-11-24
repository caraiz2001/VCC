"""
This script implements a GSEA-like algorithm for gene pairs.
For simplicity, we will not consider weights (just the rank of the pairs).
Input is a csv file with propd results containing the following columns:
- Partner: Index of the partner gene (partner > pair)
- Pair: Index of the pair gene (pair < partner)
- theta: Differential proportionality metric. Between 0 and 1. 
- Other columns are ignored. 

The idea is to compute a genewise 'differential expression - like' score 
based on the pairs that each gene is involved in.
"""


# ----------------------------
# Import libraries
# ----------------------------
import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
import numpy as np


# ----------------------------
# Functions
# ----------------------------

# Sorting function for big propd csv files
def sort_csv_by_score(input_csv: Path, output_sorted_csv: Path):
    """
    Produces a sorted CSV/TSV with the SAME columns ("Partner","Pair","theta","FDR"),
    ordered by the theta column (3rd column), ascending (close to 0 first).
    Detects delimiter (comma or tab) automatically.
    Requires coreutils: tail, head, sort.
    """
    # Detect delimiter from header
    with open(input_csv, "r") as f:
        header = f.readline()
        delimiter = "\t" if "\t" in header else ","
        sort_flag = "-t$'\\t'" if delimiter == "\t" else "-t,"

    tmpdir = tempfile.mkdtemp(prefix="pairgsea_")
    tmp_nohdr = Path(tmpdir) / "nohdr.txt"
    tmp_sorted = Path(tmpdir) / "sorted.txt"

    # Strip header
    with open(tmp_nohdr, "wb") as out:
        subprocess.run(["tail", "-n", "+2", str(input_csv)], stdout=out, check=True)

    # Numeric sort on 3rd column (theta), ascending
    sort_key = "-k3,3g"
    sort_cmd = ["sort", sort_flag, sort_key, str(tmp_nohdr)]
    with open(tmp_sorted, "wb") as out:
        subprocess.run(" ".join(sort_cmd), shell=True, stdout=out, check=True)

    # Reattach header
    with open(output_sorted_csv, "wb") as out:
        subprocess.run(["head", "-n", "1", str(input_csv)], stdout=out, check=True)
        subprocess.run(["cat", str(tmp_sorted)], stdout=out, check=True)

# Function to compute enrichment score for each gene
def compute_es_stream(sorted_csv: Path,
                      gene_id: dict, N: int):
    """
    Single streaming pass over the sorted file (rank order).
    Updates each gene's running-sum:
      - misses: constant decrement = 1/(N - k_g)
      - hits:   +1/k_g (unweighted) OR +w_i / sum_w_g (weighted)
    Tracks max and min to get ES+ and ES-.
    """
    G = len(gene_id) # G is the number of genes
    K = G - 1 # K is the 'set size', which equals the number of pairs per gene = G - 1

    if N == (G * (G - 1)) / 2:
        print(f"Input file looks correct: {N} pairs for {G} genes.")
    else:
        raise RuntimeError(f"Input file has {N} pairs, but expected {G * (G - 1) / 2} for {G} genes.")

    # Precompute increments for misses and hits.
    # As the gene set size is constant, the miss is constant too
    miss_inc = 1.0 / (N - K) # Pmiss = 1 / (total pairs - pairs with gene g) = 1 / [(G*G - G) - (G-1)]
    hit_const = 1.0 / K # Phit = 1 / (pairs with gene g) = 1 / (G-1)
    
    # States
    last_idx = np.zeros(G, dtype=np.int64)  # last position seen; start at 0
    cur = np.zeros(G, dtype=np.float64) # current running sum
    maxv = np.zeros(G, dtype=np.float64) # max ES value
    minv = np.zeros(G, dtype=np.float64) # min ES value (negative)

    # Detect delimiter from header
    with open(sorted_csv, newline="") as f:
        header = f.readline()
        delimiter = "\t" if "\t" in header else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delimiter)
        i = 0
        for row in rdr:
            i += 1
            g1 = int(row["Partner"].strip())
            g2 = int(row["Pair"].strip())
            # genes must exist from prepass
            id1 = gene_id[g1]
            id2 = gene_id[g2]

            # Unweighted fast path
            dmiss = (i - last_idx[id1] - 1) * miss_inc  # miss_inc is scalar
            if dmiss:
                cur[id1] -= dmiss
                if cur[id1] < minv[id1]: minv[id1] = cur[id1]
            cur[id1] += hit_const
            if cur[id1] > maxv[id1]: maxv[id1] = cur[id1]
            last_idx[id1] = i
            # gene 2
            dmiss = (i - last_idx[id2] - 1) * miss_inc
            if dmiss:
                cur[id2] -= dmiss
                if cur[id2] < minv[id2]: minv[id2] = cur[id2]
            cur[id2] += hit_const
            if cur[id2] > maxv[id2]: maxv[id2] = cur[id2]
            last_idx[id2] = i

    # Tail drift (after last hit to end of list)
    tail = N - last_idx
    cur -= tail * miss_inc
    # update minima where tail pushes below existing min
    minv = np.minimum(minv, cur)

    # ES selection (direction with larger |ES|)
    pos = maxv
    neg = minv  # negative numbers (or zero)
    es = np.where(np.abs(pos) >= np.abs(neg), pos, neg)
    return es, pos, neg


def write_results(out_csv: Path, gene_id, es, pos, neg, pvals, qvals):
    """
    Write out results to CSV. The output file has the following columns:
    gene, ES, ES_pos, ES_neg, set_size, pval, qval
    """
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "ES", "ES_pos", "ES_neg", "set_size", "pval", "qval"])
        for i, g in enumerate(gene_id):
            w.writerow([g, f"{es[i]:.8g}", f"{pos[i]:.8g}", f"{neg[i]:.8g}",
                        f"{pvals[i]:.6g}", f"{qvals[i]:.6g}"])
def prepass_counts(input_csv: Path):
    """
    Reads CSV once (unsorted), builds:
      - gene -> id mapping (Partner and Pair are numeric IDs)
      - k[g]: number of pairs containing gene g (always G-1)
      - N: total pairs (rows)
    Requires columns "Partner","Pair","theta" (header present).
    """
    gene_ids = set()
    N = 0

    # Detect delimiter from header
    with open(input_csv, newline="") as f:
        header = f.readline()
        delimiter = "\t" if "\t" in header else ","
        f.seek(0)
        rdr = csv.DictReader(f, delimiter=delimiter)
        for row in rdr:
            g1 = int(row["Partner"].strip())
            g2 = int(row["Pair"].strip())
            gene_ids.add(g1)
            gene_ids.add(g2)
            N += 1

    if N == 0:
        raise RuntimeError("Input appears empty.")

    gene_ids = sorted(gene_ids)
    gene_id = {g: i for i, g in enumerate(gene_ids)}
    # id_to_gene = gene_ids  # list of integer gene indexes
    G = len(gene_ids)
    
    if N != (G * (G - 1)) / 2:
        raise RuntimeError(f"Input has {N} pairs, but expected {G * (G - 1) / 2} for {G} genes.")

    return gene_id, N

    return gene_id, N

def ks_one_sample_pvalue_asymptotic(D: float, n: int) -> float:
    """
    Asymptotic p-value for one-sample KS test:
    P(D_n >= D) ≈ 2 * sum_{j=1..∞} (-1)^{j-1} exp(-2 j^2 λ^2),
    with λ = D * (sqrt(n) + 0.12 + 0.11 / sqrt(n)).
    Good when n is reasonably large (true here).
    """
    if n <= 0 or D <= 0:
        return 1.0
    sqrtn = math.sqrt(n)
    lam = (sqrtn + 0.12 + 0.11 / sqrtn) * D
    # sum until terms are tiny
    total = 0.0
    j = 1
    while True:
        term = math.exp(-2.0 * (j * j) * (lam * lam))
        add = (1 if j % 2 == 1 else -1) * term
        total += add
        if term < 1e-12:  # tighten if you like
            break
        j += 1
        if j > 1000:  # absolute cap
            break
    p = 2.0 * total
    # numerical guard
    if p < 0.0:
        p = 0.0
    if p > 1.0:
        p = 1.0
    return p

# Benjamini-Hochberg FDR for multiple testing
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR (q-values). Vectorized; stable for ties.
    """
    m = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * m / (np.arange(m) + 1)
    # enforce monotonicity
    for i in range(m - 2, -1, -1):
        if q[i] > q[i + 1]:
            q[i] = q[i + 1]
    out = np.empty_like(q)
    out[order] = q
    return out




# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Pair-wise GSEA (KS enrichment per gene over ranked pairs)."
    )
    ap.add_argument("input_csv", type=Path, help="CSV with Partner,Pair,theta,FDR")
    ap.add_argument("--presorted", action="store_true",
                    help="Input is already sorted globally by the chosen score/order.")
    ap.add_argument("--sorted-out", type=Path, default=None,
                    help="Optional path to write the sorted CSV (or reuse if presorted).")
    ap.add_argument("--out", type=Path, default=Path("pair_gsea_results.csv"),
                    help="Output CSV with ES/p-values/FDR.")
    args = ap.parse_args()

    # 1) Pre-pass on unsorted file (fast, streaming)
    print("Pre-pass: counting per-gene set sizes...", file=sys.stderr)
    gene_id, N = prepass_counts(args.input_csv)

    G = len(gene_id)
    print(f"Found {G} genes; {N} pairs (rows).", file=sys.stderr)

    # 2) Sort input file by score (theta), ascending
    #    Either confirm presorted, or produce a sorted copy cheaply with GNU sort
    if args.presorted:
        sorted_csv = args.input_csv if args.sorted_out is None else args.sorted_out
        if sorted_csv != args.input_csv:
            # copy (optional)
            subprocess.run(["cp", str(args.input_csv), str(sorted_csv)], check=True)
    else:
        sorted_csv = args.sorted_out or Path("pairs.sorted.csv")
        sort_csv_by_score(args.input_csv, sorted_csv)
        print("Sorting done.", file=sys.stderr)

    # 3) Streaming ES over the ranked file
    print("Computing enrichment scores (streaming pass)...", file=sys.stderr)
    es, pos, neg = compute_es_stream(sorted_csv, gene_id, N)

    # 4) p-values (one-sample KS) and BH-FDR
    print("Computing p-values and FDR...", file=sys.stderr)
    # KS uses D = |ES| and sample size = k (hits for that gene)
    D = np.abs(es)
    K = G - 1  # set size (pairs per gene)
    pvals = np.array([ks_one_sample_pvalue_asymptotic(D[i], K) for i in range(G)], dtype=np.float64)
    qvals = bh_fdr(pvals)

    # 5) Write results
    write_results(args.out, gene_id, es, pos, neg, pvals, qvals)
    print(f"Done. Results -> {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()


