#!/usr/bin/env python3
"""
CCDS Codon Usage Explorer – data preprocessing
Reads:  data/CCDS_codon_usage.tsv
Writes: docs/data.json   – gene metadata, PCA coords, cluster labels, centroids
        docs/codons.json – raw codon-count matrix for per-gene detail views
"""
import json
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# ── Config ────────────────────────────────────────────────────────────────────
INPUT   = "data/CCDS_codon_usage.tsv"
OUT_DIR = "docs"
N_CLUSTERS = 12   # change if you want more / fewer
RANDOM_STATE = 42
# ─────────────────────────────────────────────────────────────────────────────

def log(msg):
    print(msg, flush=True)

log(f"Reading {INPUT}…")
df = pd.read_csv(INPUT, sep="\t", dtype={"chromosome": str})

# Drop trailing unnamed column produced by trailing tab
df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

# Identify codon columns
codon_cols  = [c for c in df.columns if "-" in c]
stop_codons = [c for c in codon_cols if "stop" in c.lower()]
sense_codons= [c for c in codon_cols if "stop" not in c.lower()]

log(f"  {len(df):,} genes | {len(codon_cols)} codon cols "
    f"({len(sense_codons)} sense, {len(stop_codons)} stop)")

# Fill NaN with 0 (missing entry = codon not observed)
df[codon_cols] = df[codon_cols].fillna(0).astype(float)

# Normalise: frequency = count / total sense codons (exclude stop from denominator)
sense_totals = df[sense_codons].sum(axis=1).replace(0, 1)
sense_freq   = df[sense_codons].div(sense_totals, axis=0)
codon_freq   = df[codon_cols].div(sense_totals, axis=0)   # all 64, for display

# ── PCA ───────────────────────────────────────────────────────────────────────
log("Running PCA on sense-codon frequencies…")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(sense_freq)

pca = PCA(n_components=3, random_state=RANDOM_STATE)
pcs = pca.fit_transform(X_scaled)
ev  = pca.explained_variance_ratio_
log(f"  Explained variance: PC1={ev[0]:.2%}  PC2={ev[1]:.2%}  PC3={ev[2]:.2%}")

# ── Clustering ────────────────────────────────────────────────────────────────
log(f"Running KMeans (k={N_CLUSTERS})…")
km = KMeans(n_clusters=N_CLUSTERS, random_state=RANDOM_STATE,
            n_init=10, max_iter=300)
labels = km.fit_predict(X_scaled)

cluster_sizes = {int(k): int(v)
                 for k, v in zip(*np.unique(labels, return_counts=True))}
log(f"  Cluster sizes: {cluster_sizes}")

# Centroid in frequency space (inverse-transform, clip negatives)
centroid_scaled = km.cluster_centers_
centroid_freq   = scaler.inverse_transform(centroid_scaled)
centroid_freq   = np.clip(centroid_freq, 0, None)

# Cluster centroids for all 64 codons (for heatmap – include stop codons too)
# We'll rebuild centroid for all codons from per-gene averages
cluster_centroids_all = []
for j in range(N_CLUSTERS):
    mask = labels == j
    mean_freq = codon_freq[mask].mean(axis=0).values
    cluster_centroids_all.append([round(float(v), 4) for v in mean_freq])

# ── Gene loading-bar amino-acid grouping for heatmap x-axis ──────────────────
aa_map = {}
for c in codon_cols:
    parts = c.split("-")
    aa    = parts[1] if len(parts) > 1 else "?"
    aa_map.setdefault(aa, []).append(c)

# Sort codons by amino acid for a prettier heatmap
codon_order = []
for aa in sorted(aa_map):
    codon_order.extend(sorted(aa_map[aa]))

codon_order_idx = [codon_cols.index(c) for c in codon_order]

# Reorder centroid columns
cluster_centroids_reordered = []
for row in cluster_centroids_all:
    cluster_centroids_reordered.append([row[i] for i in codon_order_idx])

# ── Build data.json ───────────────────────────────────────────────────────────
log("Building data.json…")

# Compact gene array: [gene, ccds_id, chr, length, strand, cluster, pc1, pc2, pc3]
genes = []
for i in range(len(df)):
    row = df.iloc[i]
    genes.append([
        str(row.get("gene", "")),
        str(row.get("ccds_id", "")),
        str(row.get("chromosome", "")),
        int(row.get("Length", 0)) if pd.notna(row.get("Length")) else 0,
        str(row.get("cds_strand", "")),
        int(labels[i]),
        round(float(pcs[i, 0]), 3),
        round(float(pcs[i, 1]), 3),
        round(float(pcs[i, 2]), 3),
    ])

data_out = {
    "columns":           ["gene","ccds_id","chr","len","strand","cluster","pc1","pc2","pc3"],
    "genes":             genes,
    "codon_cols":        codon_cols,      # original order (for codons.json lookup)
    "codon_cols_sorted": codon_order,     # sorted by AA (for heatmap)
    "sense_codons":      sense_codons,
    "stop_codons":       stop_codons,
    "n_clusters":        N_CLUSTERS,
    "cluster_sizes":     [cluster_sizes.get(j, 0) for j in range(N_CLUSTERS)],
    "centroids":         cluster_centroids_reordered,   # [cluster][codon_sorted]
    "explained_variance":[round(float(v), 4) for v in ev],
}

os.makedirs(OUT_DIR, exist_ok=True)
data_path = os.path.join(OUT_DIR, "data.json")
with open(data_path, "w") as f:
    json.dump(data_out, f, separators=(",", ":"))

log(f"  Wrote {data_path} ({os.path.getsize(data_path)/1e6:.1f} MB)")

# ── Build codons.json – raw count matrix ──────────────────────────────────────
log("Building codons.json (raw counts, for per-gene detail)…")

# Store as integer counts to keep file small
counts_matrix = df[codon_cols].astype(int).values.tolist()

codons_out = {"matrix": counts_matrix}
codons_path = os.path.join(OUT_DIR, "codons.json")
with open(codons_path, "w") as f:
    json.dump(codons_out, f, separators=(",", ":"))

log(f"  Wrote {codons_path} ({os.path.getsize(codons_path)/1e6:.1f} MB)")
log("Done ✓")
