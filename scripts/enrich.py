# https://string-db.org/cgi/help?subpage=api%23valuesranks-enrichment-api
import argparse
import json
import random
from pathlib import Path

import requests  # python -m pip install requests

string_api_url = "https://version-12-0.string-db.org/api"
output_format = "json"
method = "geneset_description"
request_url = "/".join([string_api_url, output_format, method])


def load_data(path):
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def group_genes_by_cluster(data):
    columns = data.get("columns", [])
    if "gene" not in columns or "cluster" not in columns:
        raise KeyError("Expected 'gene' and 'cluster' columns in metadata JSON")

    gene_idx = columns.index("gene")
    cluster_idx = columns.index("cluster")
    clusters = {}
    for row in data.get("genes", []):
        cluster = str(row[cluster_idx])
        gene = str(row[gene_idx])
        clusters.setdefault(cluster, []).append(gene)

    # Preserve input order while removing duplicates per cluster
    return {k: list(dict.fromkeys(v)) for k, v in clusters.items()}


def load_cluster_genes(path, cluster_id=None):
    data = load_data(path)
    if "clusters" in data:
        clusters = {
            str(k): list(dict.fromkeys(v))
            for k, v in data["clusters"].items()
        }
    elif "columns" in data and "genes" in data:
        clusters = group_genes_by_cluster(data)
    else:
        raise KeyError(
            "Input JSON must contain either 'clusters' or gene metadata with 'columns' and 'genes'."
        )

    if cluster_id is None:
        return clusters

    cluster_key = str(cluster_id)
    if cluster_key not in clusters:
        raise KeyError(f"Cluster {cluster_id} not found in {path}")

    return {cluster_key: clusters[cluster_key]}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run STRING gene set description on genes loaded from a cluster metadata file."
    )
    parser.add_argument(
        "--input",
        default="docs/data.json",
        help="Path to JSON file containing either cluster metadata or cluster->gene mapping.",
    )
    parser.add_argument(
        "--cluster",
        type=int,
        help="Cluster ID to use from the input file.",
    )
    parser.add_argument(
        "--species",
        type=int,
        default=9606,
        help="NCBI taxon identifier (default: 9606 for human).",
    )
    parser.add_argument(
        "--caller-identity",
        default="www.awesome_app.org",
        help="Caller identity for STRING API.",
    )
    parser.add_argument(
        "--output-file",
        default="docs/cluster_descriptions.json",
        help="Write top descriptions per cluster to this JSON file.",
    )
    return parser.parse_args()


def describe_row(row):
    return {
        "primary": row.get("primary_description", "-"),
        "secondary": row.get("secondary_description", "-"),
        "tertiary": row.get("tertiary_description", "-"),
    }


def fetch_cluster_descriptions(genes, species, caller_identity):
    identifiers = "%0d".join(genes)
    params = {
        "identifiers": identifiers,
        "species": species,
        "caller_identity": caller_identity,
    }

    response = requests.post(request_url, data=params)
    response.raise_for_status()
    data = response.json()
    return [describe_row(row) for row in data]


def main():
    args = parse_args()
    clusters = load_cluster_genes(args.input, args.cluster)
    if not clusters:
        raise SystemExit("No genes found for the requested cluster(s).")

    results = {}
    for cluster_id, genes in clusters.items():
        if not genes:
            continue

        if len(genes) > 2000:
            sampled_genes = random.Random(42).sample(genes, 2000)
        else:
            sampled_genes = genes

        print(f"Cluster {cluster_id}: original={len(genes)}, sample={len(sampled_genes)}")
        descriptions = fetch_cluster_descriptions(
            sampled_genes, args.species, args.caller_identity
        )
        top_three = descriptions[:3]
        results[cluster_id] = {
            "sampled_size": len(sampled_genes),
            "top_descriptions": top_three,
        }

        cluster_output_path = Path(args.output_file).parent / f"cluster_{cluster_id}_descriptions.json"
        cluster_output_path.parent.mkdir(parents=True, exist_ok=True)
        with cluster_output_path.open("w", encoding="utf-8") as f:
            json.dump(
                {
                    "cluster": cluster_id,
                    "sampled_size": len(sampled_genes),
                    "top_descriptions": top_three,
                },
                f,
                indent=2,
            )

        for i, desc in enumerate(top_three, start=1):
            print(f"  Primary category {i}: {desc['primary']}")
            print(f"    Secondary: {desc['secondary']}")
            print(f"    Tertiary:  {desc['tertiary']}")
        print(f"  Wrote per-cluster output to {cluster_output_path}")
        print()

    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump({"clusters": results}, f, indent=2)

    print(f"Wrote overall descriptions to {output_path}")


if __name__ == "__main__":
    main()
