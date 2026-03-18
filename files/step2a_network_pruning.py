"""
Step 2a. Map Entrez GeneIDs to symbols and prune the STRING network
to only the genes that survived the TPM filter.

Loads
-----
.cache/step1_filtered_genes.parquet

Saves
-----
.cache/step2a_filtered_genes.parquet   (with string_id column added)
.cache/step2a_graph.pkl                (full pruned NetworkX graph)
"""

import mygene
import networkx as nx
import pandas as pd
from tqdm import tqdm

import cache


# ---------------------------------------------------------------------------
# 2.1  Map Entrez GeneIDs → gene symbols
# ---------------------------------------------------------------------------

def map_entrez_to_symbols(filtered_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Query MyGene.info to add a 'symbol' column to `filtered_genes`.
    """
    print("--- Mapping Entrez IDs to Symbols ---")
    filtered_genes = filtered_genes.drop_duplicates(subset="GeneID").copy()

    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        filtered_genes["GeneID"].tolist(),
        scopes="entrezgene",
        fields="symbol",
        species="human",
        verbose=False,
    )

    mapping = {}
    for r in tqdm(results, desc="Processing mapping results"):
        if "symbol" in r:
            mapping[str(r["query"])] = r["symbol"]

    filtered_genes["symbol"] = filtered_genes["GeneID"].astype(str).map(mapping)
    print(f"✅ Symbols resolved: {filtered_genes['symbol'].notna().sum()} / {len(filtered_genes)}")
    return filtered_genes


# ---------------------------------------------------------------------------
# 2.2  Prune the STRING network to filtered genes only
# ---------------------------------------------------------------------------

def build_pruned_network(
    filtered_genes: pd.DataFrame,
    links_path: str,
    info_path: str,
) -> tuple[pd.DataFrame, nx.Graph]:
    """
    1. Merge filtered genes with STRING protein info to get string_id.
    2. Keep only STRING edges where both endpoints are in filtered genes.
    """
    print(f"--- Loading STRING data from {links_path} ---")
    links = pd.read_csv(links_path, sep=" ")
    # links = links[links["combined_score"] > 700]
    info  = pd.read_csv(info_path, sep="\t")[["#string_protein_id", "preferred_name"]]
    info.columns = ["string_id", "symbol"]

    print("Merging filtered genes with STRING info...")
    filtered_genes = filtered_genes.merge(info, on="symbol", how="inner")

    valid_ids = set(filtered_genes["string_id"])
    print(f"Pruning edges for {len(valid_ids)} valid proteins...")
    mask         = links["protein1"].isin(valid_ids) & links["protein2"].isin(valid_ids)
    pruned_links = links[mask].copy()

    print("Building NetworkX graph object...")
    G = nx.from_pandas_edgelist(
        pruned_links,
        source="protein1",
        target="protein2",
        edge_attr="combined_score",
    )
    G.add_nodes_from(valid_ids)  # include isolated nodes

    print(f"✅ Network created — Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    return filtered_genes, G


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("🚀 Starting Step 2a: Network Pruning")

    filtered_genes = cache.load_df("step1_filtered_genes.parquet")
    print(f"Loaded {len(filtered_genes)} genes from Step 1.")

    filtered_genes = map_entrez_to_symbols(filtered_genes)
    filtered_genes, G = build_pruned_network(
        filtered_genes,
        links_path="other_dataset/9606.protein.links.v12.0.txt",
        info_path="other_dataset/9606.protein.info.v12.0.txt",
    )

    print("--- Saving results to cache ---")
    cache.save_df(filtered_genes, "step2a_filtered_genes.parquet")
    cache.save_graph(G, "step2a_graph.pkl")

    print("\n--- Preview of Results ---")
    print(filtered_genes[["symbol", "string_id"]].head(10))
    print("\n🎉 Step 2a Complete.")
