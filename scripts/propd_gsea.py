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
# Helpers
# ----------------------------

# Theta should be between 0 and 1, so I don't think we need this function 
# def score_from_theta(theta: float, transform: str) -> float:
#     if transform == "raw":
#         return theta
#     elif transform == "abs":
#         return abs(theta)
#     elif transform == "absdiff1":
#         return abs(theta - 1.0)
#     else:
#         raise ValueError(f"Unknown score_transform: {transform}")

# Function to obtain an analytical pvalue, no permutations needed (see more satistical details)
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
# Pass 1: Map genes and count occurrences (and optionally weight sums)
# ----------------------------

def prepass_counts(input_csv: Path, score_transform: str, weighted_power: float):
    """
    Reads CSV once (unsorted), builds:
      - gene -> id mapping (here, Partner and Pair are already numeric IDs)
      - k[g]: number of pairs containing gene g
      - sum_w[g]: sum of weights (|score|^p) of rows where g appears (if p>0)
      - N: total pairs (rows)
    Requires columns "Partner","Pair","theta" (header present).
    Streaming, low memory.
    """
    # Since Partner and Pair are already numeric IDs, use them directly
    k = defaultdict(int)
    sum_w = defaultdict(float) if weighted_power > 0.0 else None
    N = 0
    gene_ids = set()

    with open(input_csv, newline="") as f:
        rdr = csv.DictReader(f)
        # If weighted_power=0, sum_w is equal to k (number of pairs)
        for row in rdr:
            # extracts gene indexes and theta for the row
            g1 = int(row["Partner"])
            g2 = int(row["Pair"])
            th = float(row["theta"])

            # Adds the gene indexes (gene_ids is a set, so no duplicates)
            gene_ids.add(g1)
            gene_ids.add(g2)

            # Number of pairs containing gene g. Actually, this is not necessary as the pairs are all vs all
            k[g1] += 1
            k[g2] += 1

            if sum_w is not None:
                #s = score_from_theta(th, score_transform)
                s = th

                # weighted power refers to 'p' in the paper. If it is set to 0, the ranking metric is not weighted for the score (all pairs contribute equally)
                # no need of absolute because theta is between 0 and 1
                # w = abs(s)**weighted_power

                # w = 1 if weighted_power = 0
                w = s ** weighted_power
                # Sum-w for gene A is the sum of weights of all pairs involving A (or just the number of genes if weighted_power=0)
                sum_w[g1] += w
                sum_w[g2] += w
            N += 1
            # gene_ids is already a set of integer gene indexes
            gene_ids = sorted(gene_ids)
            gene_id = {g: i for i, g in enumerate(gene_ids)}
            id_to_gene = gene_ids  # list of integer gene indexes
            G = len(gene_ids)
            k_arr = np.zeros(G, dtype=np.int64)
            for g, i in gene_id.items():
                k_arr[i] = k[g]
            if sum_w is not None:
                sumw_arr = np.zeros(G, dtype=np.float64)
                for g, i in gene_id.items():
                    val = sum_w.get(g, 0.0)
                    sumw_arr[i] = val if val > 0 else 1e-12
            else:
                sumw_arr = None

            if N == 0:
                raise RuntimeError("Input appears empty.")

    return gene_id, id_to_gene, k_arr, sumw_arr, N

# I think all of this could be simplified because 
# K is just the number of pairs containing gene g, which is the number of genes - 1 (all vs all)
# sum_w is just K if weighted_power = 0 (also the number of genes - 1)

# ----------------------------
# Sorting: use GNU sort (fast, external, low RAM)
# ----------------------------

def sort_csv_by_score(input_csv: Path, output_sorted_csv: Path, order: str):
    """
    Produces a sorted CSV with the SAME columns ("Partner","Pair","theta","FDR"),
    ordered by the theta column (3rd column).
    Requires coreutils: tail, head, sort.
    """
    assert order in ("asc", "desc")
    tmpdir = tempfile.mkdtemp(prefix="pairgsea_")
    tmp_nohdr = Path(tmpdir) / "nohdr.csv"
    tmp_sorted = Path(tmpdir) / "sorted.csv"

    # Strip header
    with open(tmp_nohdr, "wb") as out:
        subprocess.run(["tail", "-n", "+2", str(input_csv)], stdout=out, check=True)

    # Numeric sort on 3rd column (theta)
    sort_key = f"-k3,3g{'r' if order=='desc' else ''}"
    with open(tmp_sorted, "wb") as out:
        subprocess.run(["sort", "-t,", sort_key, str(tmp_nohdr)], stdout=out, check=True)

    # Reattach header
    with open(output_sorted_csv, "wb") as out:
        subprocess.run(["head", "-n", "1", str(input_csv)], stdout=out, check=True)
        subprocess.run(["cat", str(tmp_sorted)], stdout=out, check=True)

def compute_es_stream(sorted_csv: Path,
                      gene_id: dict,
                      k_arr: np.ndarray,
                      sumw_arr: np.ndarray,
                      N: int,
                      score_transform: str,
                      weighted_power: float):
    """
    Single streaming pass over the sorted file (rank order).
    Updates each gene's running-sum:
      - misses: constant decrement = 1/(N - k_g)
      - hits:   +1/k_g (unweighted) OR +w_i / sum_w_g (weighted)
    Tracks max and min to get ES+ and ES-.
    """
    G = k_arr.shape[0] # G is the number of genes
    K = G - 1
    
    # Precompute increments
    # miss_inc = 1.0 / (N - k_arr.astype(np.float64)) 
    miss_inc = 1.0 / (N - K) # Pmiss = 1 / (total pairs - pairs with gene g) = 1 / [(G*G - G) - (G-1)]
    hit_const = np.zeros(G, dtype=np.float64)
    
    # For now, this will be False because weighted_power = 0
    weighted = weighted_power > 0.0
    
    # If not weighted, hit increment is constant = 1/k for each gene
    if not weighted:
        # hit_const = 1.0 / k_arr.astype(np.float64)
        hit_const = 1.0 / K
    else:
        if sumw_arr is None:
            raise RuntimeError("Weighted mode requires sum of weights per gene.")
    
    # States
    last_idx = np.zeros(G, dtype=np.int64) # last position seen; start at 0
    cur = np.zeros(G, dtype=np.float64)
    maxv = np.zeros(G, dtype=np.float64)
    minv = np.zeros(G, dtype=np.float64)

    # Streaming over ranked rows
    i = 0
    with open(sorted_csv, newline="") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            i += 1
            g1 = row["Partner"]
            g2 = row["Pair"]
            # genes must exist from prepass
            id1 = gene_id[g1]
            id2 = gene_id[g2]

            # For now weighted is False 
            if weighted:
                th = float(row["theta"])
                #s = score_from_theta(th, score_transform)
                s = th
                w = abs(s)**weighted_power
                # Update gene 1
                # apply drift from intervening misses
                
                dmiss = (i - last_idx[id1] - 1) * miss_inc[id1] 
                if dmiss:
                    cur[id1] -= dmiss
                    if cur[id1] < minv[id1]: minv[id1] = cur[id1]
                cur[id1] += (w / sumw_arr[id1])
                if cur[id1] > maxv[id1]: maxv[id1] = cur[id1]
                last_idx[id1] = i

                # Update gene 2
                dmiss = (i - last_idx[id2] - 1) * miss_inc[id2]
                if dmiss:
                    cur[id2] -= dmiss
                    if cur[id2] < minv[id2]: minv[id2] = cur[id2]
                cur[id2] += (w / sumw_arr[id2])
                if cur[id2] > maxv[id2]: maxv[id2] = cur[id2]
                last_idx[id2] = i
            else:
                # Unweighted fast path
                # Calculate the miss penalty
                # For gene in Partner
                dmiss = (i - last_idx[id1] - 1) * miss_inc[id1] # If the last index is the previous row, dmiss = 0
                if dmiss:
                    cur[id1] -= dmiss
                    if cur[id1] < minv[id1]: minv[id1] = cur[id1]
                cur[id1] += hit_const[id1]
                if cur[id1] > maxv[id1]: maxv[id1] = cur[id1]
                last_idx[id1] = i
                # gene 2
                dmiss = (i - last_idx[id2] - 1) * miss_inc[id2]
                if dmiss:
                    cur[id2] -= dmiss
                    if cur[id2] < minv[id2]: minv[id2] = cur[id2]
                cur[id2] += hit_const[id2]
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

def write_results(out_csv: Path, id_to_gene, es, pos, neg, k_arr, pvals, qvals):
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "ES", "ES_pos", "ES_neg", "set_size", "pval", "qval"])
        for i, g in enumerate(id_to_gene):
            w.writerow([g, f"{es[i]:.8g}", f"{pos[i]:.8g}", f"{neg[i]:.8g}",
                        int(k_arr[i]), f"{pvals[i]:.6g}", f"{qvals[i]:.6g}"])

# ----------------------------
# CLI
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Pair-wise GSEA (KS enrichment per gene over ranked pairs)."
    )
    ap.add_argument("input_csv", type=Path, help="CSV with Partner,Pair,theta,FDR")
    ap.add_argument("--score-transform", choices=["raw","abs","absdiff1"], default="raw",
                    help="Ranking score from theta (default: raw).")
    ap.add_argument("--order", choices=["asc","desc"], default="asc",
                    help="Sort order for the score (default: asc).")
    ap.add_argument("--weighted-power", type=float, default=0.0,
                    help="Use |score|^p weights on hits (p=0 => unweighted).")
    ap.add_argument("--presorted", action="store_true",
                    help="Input is already sorted globally by the chosen score/order.")
    ap.add_argument("--sorted-out", type=Path, default=None,
                    help="Optional path to write the sorted CSV (or reuse if presorted).")
    ap.add_argument("--out", type=Path, default=Path("pair_gsea_results.csv"),
                    help="Output CSV with ES/p-values/FDR.")
    args = ap.parse_args()

    # 1) Pre-pass on unsorted file (fast, streaming)
    print("Pre-pass: counting per-gene set sizes...", file=sys.stderr)
    gene_id, id_to_gene, k_arr, sumw_arr, N = prepass_counts(
        args.input_csv, args.score_transform, args.weighted_power
    )
    G = len(id_to_gene)
    print(f"Found {G} genes; {N} pairs (rows).", file=sys.stderr)

    # 2) Sorting: either confirm presorted, or produce a sorted copy cheaply with GNU sort
    if args.presorted:
        sorted_csv = args.input_csv if args.sorted_out is None else args.sorted_out
        if sorted_csv != args.input_csv:
            # copy (optional)
            subprocess.run(["cp", str(args.input_csv), str(sorted_csv)], check=True)
    else:
        sorted_csv = args.sorted_out or Path("pairs.sorted.csv")
        print(f"Sorting by {args.score_transform} ({args.order}) into {sorted_csv} ...", file=sys.stderr)
        sort_csv_by_score(args.input_csv, sorted_csv, args.score_transform, args.order)
        print("Sorting done.", file=sys.stderr)

    # 3) Streaming ES over the ranked file
    print("Computing enrichment scores (streaming pass)...", file=sys.stderr)
    es, pos, neg = compute_es_stream(sorted_csv, gene_id, k_arr, sumw_arr, N,
                                     args.score_transform, args.weighted_power)

    # 4) p-values (one-sample KS) and BH-FDR
    print("Computing p-values and FDR...", file=sys.stderr)
    # KS uses D = |ES| and sample size = k (hits for that gene)
    D = np.abs(es)
    pvals = np.array([ks_one_sample_pvalue_asymptotic(D[i], int(k_arr[i])) for i in range(G)], dtype=np.float64)
    qvals = bh_fdr(pvals)

    # 5) Write results
    write_results(args.out, id_to_gene, es, pos, neg, k_arr, pvals, qvals)
    print(f"Done. Results -> {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
