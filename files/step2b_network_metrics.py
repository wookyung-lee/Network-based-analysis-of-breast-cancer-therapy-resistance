"""
Step 2b. Compute Node Degree (ND) and Betweenness Centrality (BC)
for all nodes in the pruned network.

Definitions
-----------
Node Degree            : number of edges connected to a node.
Betweenness Centrality : fraction of shortest paths (between all node pairs)
                         that pass through a given node.

Loads
-----
.cache/step2a_filtered_genes.parquet
.cache/step2a_graph.pkl

Saves
-----
.cache/step2b_filtered_genes.parquet   (with ND and BC columns added)
"""

import networkx as nx
import pandas as pd

import cache


# ---------------------------------------------------------------------------
# 2.3  Compute ND and BC
# ---------------------------------------------------------------------------

def compute_network_metrics(
    filtered_genes: pd.DataFrame,
    G: nx.Graph,
    bc_k: int = 500,
) -> pd.DataFrame:
    """
    Attach Node Degree (ND) and Betweenness Centrality (BC) to filtered_genes.

    Parameters
    ----------
    bc_k : pivot nodes sampled for approximate BC. None = exact (slow).
    """
    print("--- Computing Network Metrics ---")

    print("Calculating Node Degree...")
    nd = dict(G.degree())

    print(f"Calculating Betweenness Centrality (k={bc_k}, this may take a minute)...")
    bc = (
        nx.betweenness_centrality(G, normalized=True)
        if bc_k is None
        else nx.betweenness_centrality(G, k=bc_k, normalized=True, seed = 42)
    )

    filtered_genes = filtered_genes.copy()
    filtered_genes["ND"] = filtered_genes["string_id"].map(nd).fillna(0)
    filtered_genes["BC"] = filtered_genes["string_id"].map(bc).fillna(0)

    print("✅ Metrics computation complete.")
    return filtered_genes


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("🚀 Starting Step 2b: Node Degree & Betweenness Centrality")

    filtered_genes = cache.load_df("step2a_filtered_genes.parquet")
    G              = cache.load_graph("step2a_graph.pkl")
    print(f"Loaded {len(filtered_genes)} genes and graph with {G.number_of_nodes()} nodes.")

    filtered_genes = compute_network_metrics(filtered_genes, G)

    print("--- Saving results to cache ---")
    cache.save_df(filtered_genes, "step2b_filtered_genes.parquet")

    print("\n--- Top 10 Genes by Betweenness Centrality ---")
    print(filtered_genes[["symbol", "ND", "BC"]].sort_values("BC", ascending=False).head(10))
    print("\n🎉 Step 2b Complete.")
