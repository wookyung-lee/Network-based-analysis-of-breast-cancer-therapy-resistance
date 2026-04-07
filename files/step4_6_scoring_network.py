"""
Steps 4 – 6.  Score genes, rank them, and build the core network.

Score function (Step 4)
-----------------------
S(n_i) = 0.50 * |log2FC(n_i)| / max(|log2FC|)
        + 0.25 * ND(n_i)     / max(ND)
        + 0.25 * BC(n_i)     / max(BC)

Purpose: rank genes by their combined biological relevance (differential
expression) and network importance (connectivity and centrality).

Step 5  — Select the top-100 nodes by score.
Step 6  — Extract the subgraph of those nodes; prune isolated singletons.
           If the connected subgraph has fewer than 50 nodes, expand the
           candidate set in increments of 25 until the threshold is met.

Loads
-----
.cache/step3_filtered_genes.parquet
.cache/step2_graph.pkl

Saves
-----
.cache/step4_6_ranked_genes.parquet
.cache/step4_6_core_net.pkl
"""

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

import cache


# ---------------------------------------------------------------------------
# Step 4 — Score
# ---------------------------------------------------------------------------

def score_genes(filtered_genes: pd.DataFrame) -> pd.DataFrame:
    """Add a 'score' column using the weighted formula above."""
    abs_fc = filtered_genes["log2FC"].abs()

    max_fc = abs_fc.max() or 1
    max_nd = filtered_genes["ND"].max() or 1
    max_bc = filtered_genes["BC"].max() or 1

    filtered_genes = filtered_genes.copy()
    filtered_genes["score"] = (
        0.50 * abs_fc                   / max_fc
        + 0.25 * filtered_genes["ND"]   / max_nd
        + 0.25 * filtered_genes["BC"]   / max_bc
    )

    print(f"Score range: {filtered_genes['score'].min():.4f} – {filtered_genes['score'].max():.4f}")
    return filtered_genes


# ---------------------------------------------------------------------------
# Step 5 — Rank and select top-N
# ---------------------------------------------------------------------------

def rank_genes(filtered_genes: pd.DataFrame) -> pd.DataFrame:
    """Return a copy sorted by score (descending) with a 'rank' column."""
    ranked = (
        filtered_genes
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    ranked["rank"] = range(1, len(ranked) + 1)
    return ranked


# ---------------------------------------------------------------------------
# Step 6 — Core network
# ---------------------------------------------------------------------------

def build_core_network(
    ranked: pd.DataFrame,
    G: nx.Graph,
    n_start: int = 100,
    min_nodes: int = 50,
    step: int = 25,
) -> tuple[nx.Graph, pd.DataFrame, int]:
    """
    Build the smallest connected subgraph with ≥ `min_nodes` nodes.
    Expands the candidate set by `step` genes at a time as needed.
    """
    n_select = n_start

    while True:
        top_genes = ranked.head(n_select).copy()
        top_ids   = set(top_genes["string_id"].dropna())

        sub        = G.subgraph(top_ids).copy()
        singletons = [n for n, deg in sub.degree() if deg == 0]
        sub.remove_nodes_from(singletons)

        n_nodes = sub.number_of_nodes()
        n_edges = sub.number_of_edges()
        print(f"n_select={n_select:4d}: connected nodes={n_nodes:4d}, edges={n_edges}")

        if n_nodes >= min_nodes:
            break
        n_select += step
        if n_select > len(ranked):
            print("Warning: exhausted all ranked genes — using full set.")
            break

    print(f"\nFinal core network: {n_nodes} nodes, {n_edges} edges")
    return sub, top_genes, n_select


# ---------------------------------------------------------------------------
# Step 6.1 — Visualise the core network
# ---------------------------------------------------------------------------

def plot_core_network(
    core_net: nx.Graph,
    id_to_symbol: dict,
    id_to_log2fc: dict,
    id_to_score: dict,
    save_path: str | None = "core_network.png",
) -> None:
    """
    Spring-layout visualisation of the core network.
    Node colour : log2FC  (blue = Resistant, red = Sensitive)
    Node size   : score   (larger = more important)
    """
    node_list   = list(core_net.nodes())
    node_colors = [id_to_log2fc.get(n, 0) for n in node_list]
    node_sizes  = [300 + 1500 * id_to_score.get(n, 0) for n in node_list]
    max_abs     = max(abs(c) for c in node_colors) or 1

    fig, ax = plt.subplots(figsize=(16, 16))
    pos = nx.spring_layout(core_net, seed=42, k=0.8)

    sc = nx.draw_networkx_nodes(
        core_net, pos,
        nodelist=node_list,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=plt.cm.coolwarm,
        vmin=-max_abs, vmax=max_abs,
        alpha=0.9, ax=ax,
    )
    nx.draw_networkx_edges(core_net, pos, alpha=0.20, width=0.7, ax=ax)
    nx.draw_networkx_labels(
        core_net, pos,
        labels={n: id_to_symbol.get(n, n) for n in node_list},
        font_size=6, font_weight="bold", ax=ax,
    )

    plt.colorbar(sc, ax=ax, label="log2FC  (red = higher in Sensitive, blue = higher in Resistant)")
    print(f"Core Network — {core_net.number_of_nodes()} nodes, {core_net.number_of_edges()} edges")
    # ax.set_title(
    #     f"Core Network — {core_net.number_of_nodes()} nodes, "
    #     f"{core_net.number_of_edges()} edges\nNode size ∝ score",
    #     fontsize=13,
    # )
    ax.axis("off")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Load step 2b and step 3 outputs
    filtered_genes = cache.load_df("step3_filtered_genes.parquet")
    G              = cache.load_graph("step2a_graph.pkl")

    filtered_genes = score_genes(filtered_genes)
    ranked         = rank_genes(filtered_genes)

    print("\nTop 10 genes:")
    print(ranked[["rank", "symbol", "log2FC", "ND", "BC", "score"]].head(10).to_string(index=False))

    core_net, _, _ = build_core_network(ranked, G)

    # Save step 4-6 outputs
    cache.save_df(ranked, "step4_6_ranked_genes.parquet")
    cache.save_graph(core_net, "step4_6_core_net.pkl")

    id_to_symbol = dict(zip(filtered_genes["string_id"], filtered_genes["symbol"]))
    id_to_score  = dict(zip(filtered_genes["string_id"], filtered_genes["score"]))
    id_to_log2fc = dict(zip(filtered_genes["string_id"], filtered_genes["log2FC"]))

    print("\nGenes in core network:", sorted(id_to_symbol.get(n, n) for n in core_net.nodes()))
    plot_core_network(core_net, id_to_symbol, id_to_log2fc, id_to_score)


# ---------------------------------------------------------------------------
# Step 6.2 — Significant-DEG network
# ---------------------------------------------------------------------------

def build_sig_network(
    sig_genes: pd.DataFrame,
    G: nx.Graph,
) -> nx.Graph:
    """
    Build a subgraph from only the significant DEGs (volcano points).
    Isolated nodes are pruned. No minimum size requirement — we take
    whatever the data gives us.

    Parameters
    ----------
    sig_genes : DataFrame from step3 get_significant_genes()
                must have 'string_id' and 'direction' columns
    G         : full pruned STRING graph from step 2a

    Returns
    -------
    sig_net : pruned NetworkX subgraph
    """
    sig_ids    = set(sig_genes["string_id"].dropna())
    sub        = G.subgraph(sig_ids).copy()
    singletons = [n for n, deg in sub.degree() if deg == 0]
    sub.remove_nodes_from(singletons)

    print(f"Significant DEG network: {sub.number_of_nodes()} nodes, "
          f"{sub.number_of_edges()} edges "
          f"({len(singletons)} singletons removed)")
    return sub


def plot_sig_network(
    sig_net: nx.Graph,
    sig_genes: pd.DataFrame,
    id_to_symbol: dict,
    save_path: str | None = "sig_network.png",
) -> None:
    """
    Visualise the significant-DEG network.

    Node colour : direction — red = up in Sensitive, blue = up in Resistant
    Node size   : |log2FC|  (larger = bigger fold change)
    """
    # Build lookup dicts from sig_genes
    id_to_fc        = dict(zip(sig_genes["string_id"], sig_genes["log2FC"]))
    id_to_direction = dict(zip(sig_genes["string_id"], sig_genes["direction"]))

    node_list  = list(sig_net.nodes())
    node_colors = [
        "red" if id_to_direction.get(n) == "up_in_sensitive" else "steelblue"
        for n in node_list
    ]
    node_sizes = [
        200 + 300 * abs(id_to_fc.get(n, 0))
        for n in node_list
    ]

    fig, ax = plt.subplots(figsize=(14, 14))
    pos = nx.spring_layout(sig_net, seed=42, k=1.0)

    nx.draw_networkx_nodes(
        sig_net, pos,
        nodelist=node_list,
        node_color=node_colors,
        node_size=node_sizes,
        alpha=0.85, ax=ax,
    )
    nx.draw_networkx_edges(sig_net, pos, alpha=0.25, width=0.8, ax=ax)
    nx.draw_networkx_labels(
        sig_net, pos,
        labels={n: id_to_symbol.get(n, n) for n in node_list},
        font_size=7, font_weight="bold", ax=ax,
    )

    # Manual legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="red",       label="Up in Sensitive"),
        Patch(facecolor="steelblue", label="Up in Resistant"),
    ]
    # ax.legend(handles=legend_elements, loc="upper left", fontsize=10)
    print(f"Significant DEG Network — {sig_net.number_of_nodes()} nodes, {sig_net.number_of_edges()} edges")
    # ax.set_title(
    #     f"Significant DEG Network — {sig_net.number_of_nodes()} nodes, "
    #     f"{sig_net.number_of_edges()} edges\nNode size ∝ |log2FC|",
    #     fontsize=13,
    # )
    ax.axis("off")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Re-run entry point to also build the sig network (append to existing run)
# ---------------------------------------------------------------------------

def run_sig_network():
    """Call this after the main scoring pipeline to add the sig-DEG network."""
    sig_genes = cache.load_df("step3_sig_genes.parquet")
    G         = cache.load_graph("step2a_graph.pkl")
    ranked    = cache.load_df("step4_6_ranked_genes.parquet")

    # Need string_id in sig_genes — merge from ranked
    if "string_id" not in sig_genes.columns:
        id_map    = ranked[["symbol", "string_id"]].drop_duplicates()
        sig_genes = sig_genes.merge(id_map, on="symbol", how="left")

    sig_net = build_sig_network(sig_genes, G)
    cache.save_graph(sig_net, "step4_6_sig_net.pkl")

    id_to_symbol = dict(zip(ranked["string_id"], ranked["symbol"]))
    plot_sig_network(sig_net, sig_genes, id_to_symbol)

    return sig_net


if __name__ == "__main__":
    # ── existing scored core network ──────────────────────────────────────────
    filtered_genes = cache.load_df("step3_filtered_genes.parquet")
    G              = cache.load_graph("step2a_graph.pkl")

    filtered_genes = score_genes(filtered_genes)
    ranked         = rank_genes(filtered_genes)

    print("\nTop 10 genes:")
    print(ranked[["rank", "symbol", "log2FC", "ND", "BC", "score"]].head(10).to_string(index=False))

    core_net, _, _ = build_core_network(ranked, G)

    cache.save_df(ranked, "step4_6_ranked_genes.parquet")
    cache.save_graph(core_net, "step4_6_core_net.pkl")

    id_to_symbol = dict(zip(filtered_genes["string_id"], filtered_genes["symbol"]))
    id_to_score  = dict(zip(filtered_genes["string_id"], filtered_genes["score"]))
    id_to_log2fc = dict(zip(filtered_genes["string_id"], filtered_genes["log2FC"]))

    print("\nGenes in core network:", sorted(id_to_symbol.get(n, n) for n in core_net.nodes()))
    plot_core_network(core_net, id_to_symbol, id_to_log2fc, id_to_score)

    # ── significant DEG network ───────────────────────────────────────────────
    print("\n--- Building Significant DEG Network ---")
    run_sig_network()
