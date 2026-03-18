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
    gene_list = filtered_genes["GeneID"].tolist()
    
    # querymany is efficient, but we can wrap the result processing
    results = mg.querymany(
        gene_list,
        scopes="entrezgene",
        fields="symbol",
        species="human",
        verbose=False
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
    1. Merge filtered genes with STRING protein info.
    2. Keep only STRING edges where both endpoints are in filtered genes.
    """
    print(f"--- Loading STRING data from {links_path} ---")
    
    # Loading large CSVs can be slow; using a simple print to confirm start
    links = pd.read_csv(links_path, sep=" ")
    info = pd.read_csv(info_path, sep="\t")[["#string_protein_id", "preferred_name"]]
    info.columns = ["string_id", "symbol"]

    print("Filtering STRING info for surviving genes...")
    filtered_genes = filtered_genes.merge(info, on="symbol", how="inner")

    valid_ids = set(filtered_genes["string_id"])
    
    print(f"Pruning edges for {len(valid_ids)} valid proteins...")
    # Filtering a large DF can be slow; no easy tqdm here without chunking, 
    # but the logic is vectorized.
    mask = links["protein1"].isin(valid_ids) & links["protein2"].isin(valid_ids)
    pruned_links = links[mask].copy()

    print("Building NetworkX graph object...")
    G = nx.from_pandas_edgelist(
        pruned_links,
        source="protein1",
        target="protein2",
        edge_attr="combined_score",
    )
    
    # Ensure all filtered genes are represented as nodes, even if isolated
    G.add_nodes_from(valid_ids)

    print(f"✅ Network created — Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    return filtered_genes, G


# ---------------------------------------------------------------------------
# 2.3  Compute ND and BC
# ---------------------------------------------------------------------------

def compute_network_metrics(
    filtered_genes: pd.DataFrame,
    G: nx.Graph,
    bc_k: int = 500,
) -> pd.DataFrame:
    """
    Attach Node Degree (ND) and Betweenness Centrality (BC).
    """
    print(f"--- Computing Network Metrics (BC sampling k={bc_k}) ---")
    
    # Degree is O(V), very fast
    print("Calculating Node Degree...")
    nd = dict(G.degree())
    
    # Betweenness Centrality is O(V*E), very slow. 
    # We use k-sampling to speed it up.
    print(f"Calculating Betweenness Centrality (this may take a minute)...")
    bc = nx.betweenness_centrality(G, k=bc_k, normalized=True)

    filtered_genes = filtered_genes.copy()
    filtered_genes["ND"] = filtered_genes["string_id"].map(nd).fillna(0)
    filtered_genes["BC"] = filtered_genes["string_id"].map(bc).fillna(0)

    print("✅ Metrics computation complete.")
    return filtered_genes


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("🚀 Starting Step 2: Protein-Interaction Network Analysis")
    
    # Load step 1 output
    filtered_genes = cache.load_df("step1_filtered_genes.parquet")
    print(f"Loaded {len(filtered_genes)} genes from Step 1.")

    filtered_genes = map_entrez_to_symbols(filtered_genes)
    
    filtered_genes, G = build_pruned_network(
        filtered_genes,
        links_path="step2/9606.protein.links.v12.0.txt",
        info_path="step2/9606.protein.info.v12.0.txt",
    )
    
    filtered_genes = compute_network_metrics(filtered_genes, G)

    # Save step 2 outputs
    print("--- Saving results to cache ---")
    cache.save_df(filtered_genes, "step2_filtered_genes.parquet")
    cache.save_graph(G, "step2_graph.pkl")

    print("\n--- Preview of Results ---")
    print(filtered_genes[["symbol", "ND", "BC"]].sort_values("BC", ascending=False).head(10))
    print("\n🎉 Step 2 Complete.")
