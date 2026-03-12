"""
Step 8. Find drugs targeting genes in the core network using DGIdb.

Queries the DGIdb GraphQL API in batches to avoid URL-size limits.
Saves all interactions and approved-drug interactions to CSV.
Generates:
  - a heatmap  : approved drugs × top-20 most-targeted genes
  - a bar chart: number of approved drugs per core-network gene

Loads
-----
.cache/step4_6_core_net.pkl
.cache/step4_6_ranked_genes.parquet   (for id → symbol mapping)

Saves
-----
drug_gene_interactions_all.csv
drug_gene_interactions_approved.csv
"""

import time

import matplotlib.pyplot as plt
import pandas as pd
import requests

import cache

DGIDB_URL = "https://dgidb.org/api/graphql"

QUERY = """
query getDrugs($genes: [String!]!) {
  genes(names: $genes) {
    nodes {
      name
      interactions {
        drug {
          name
          approved
        }
        interactionScore
      }
    }
  }
}
"""


# ---------------------------------------------------------------------------
# 8.1  Query DGIdb
# ---------------------------------------------------------------------------

def query_dgidb(
    gene_list: list[str],
    batch_size: int = 50,
    pause: float = 0.5,
) -> pd.DataFrame:
    """
    Query DGIdb for all drug-gene interactions for the supplied genes.

    Parameters
    ----------
    gene_list  : list of HGNC gene symbols
    batch_size : genes per API request (default 50)
    pause      : seconds to wait between requests (rate-limiting courtesy)

    Returns
    -------
    DataFrame with columns: gene, drug, approved, interaction_score
    """
    all_records = []
    n_batches   = (len(gene_list) + batch_size - 1) // batch_size

    for i in range(0, len(gene_list), batch_size):
        batch    = gene_list[i : i + batch_size]
        batch_no = i // batch_size + 1
        print(f"  Batch {batch_no}/{n_batches}: {len(batch)} genes ...", end=" ")

        try:
            resp = requests.post(
                DGIDB_URL,
                json={"query": QUERY, "variables": {"genes": batch}},
                timeout=60,
            )
            if resp.status_code != 200:
                print(f"HTTP {resp.status_code} — skipping")
                continue

            nodes = (
                resp.json()
                .get("data", {})
                .get("genes", {})
                .get("nodes", [])
            )

            batch_hits = 0
            for gene_node in nodes:
                for ix in gene_node.get("interactions", []):
                    drug_info = ix.get("drug") or {}
                    all_records.append({
                        "gene"             : gene_node["name"],
                        "drug"             : drug_info.get("name", ""),
                        "approved"         : drug_info.get("approved", None),
                        "interaction_score": ix.get("interactionScore", None),
                    })
                    batch_hits += 1
            print(f"{batch_hits} interactions found")

        except Exception as e:
            print(f"Error: {e} — skipping batch")

        time.sleep(pause)

    dgi_df = pd.DataFrame(all_records)
    print(f"\nTotal drug-gene pairs found: {len(dgi_df)}")
    return dgi_df


# ---------------------------------------------------------------------------
# 8.2  Summarise and plot
# ---------------------------------------------------------------------------

def summarise_drug_interactions(dgi_df: pd.DataFrame) -> pd.DataFrame:
    """Print a summary and return only approved-drug rows."""
    if dgi_df.empty:
        print("No interactions found. Verify that gene symbols are official HGNC symbols.")
        return pd.DataFrame()

    approved_df = (
        dgi_df[dgi_df["approved"] == True]
        .copy()
        .sort_values("interaction_score", ascending=False)
    )

    print(f"Total drug-gene pairs          : {len(dgi_df)}")
    print(f"Approved-drug pairs            : {len(approved_df)}")
    print(f"Unique approved drugs          : {approved_df['drug'].nunique()}")
    print(f"Core-network genes with a drug : {approved_df['gene'].nunique()}")
    print("\nTop 20 (by interaction score):")
    print(approved_df[["gene", "drug", "interaction_score"]].head(20).to_string(index=False))

    return approved_df


def plot_drug_gene_heatmap(
    approved_df: pd.DataFrame,
    top_n_genes: int = 20,
    save_path: str | None = "drug_gene_heatmap.png",
) -> None:
    """Heatmap: approved drugs × most-targeted core-network genes."""
    if approved_df.empty:
        return

    top_targeted = (
        approved_df.groupby("gene")["drug"]
        .nunique()
        .nlargest(top_n_genes)
        .index.tolist()
    )
    heat_data = (
        approved_df[approved_df["gene"].isin(top_targeted)]
        .pivot_table(index="gene", columns="drug",
                     values="interaction_score", aggfunc="max")
        .fillna(0)
    )

    fig_h = max(5,  0.5 + len(heat_data) * 0.45)
    fig_w = min(28, 1   + len(heat_data.columns) * 0.45)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat_data.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(heat_data.columns)))
    ax.set_xticklabels(heat_data.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(heat_data.index)))
    ax.set_yticklabels(heat_data.index, fontsize=9)
    plt.colorbar(im, ax=ax, label="Interaction score")
    ax.set_title(f"Approved Drugs × Core-Network Genes (top {top_n_genes} targeted genes)")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


def plot_drugs_per_gene(
    approved_df: pd.DataFrame,
    top_n: int = 25,
    save_path: str | None = "drugs_per_gene.png",
) -> None:
    """Bar chart: number of approved drugs per core-network gene."""
    if approved_df.empty:
        return

    drug_counts = (
        approved_df.groupby("gene")["drug"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_n)
    )

    fig, ax = plt.subplots(figsize=(8, max(4, len(drug_counts) * 0.38)))
    ax.barh(drug_counts.index[::-1], drug_counts.values[::-1], color="steelblue")
    ax.set_xlabel("Number of approved drugs")
    ax.set_title("Core-Network Genes by Number of Approved Drug Interactions")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


def save_drug_tables(
    dgi_df: pd.DataFrame,
    approved_df: pd.DataFrame,
    all_path: str = "drug_gene_interactions_all.csv",
    approved_path: str = "drug_gene_interactions_approved.csv",
) -> None:
    dgi_df.to_csv(all_path, index=False)
    approved_df.to_csv(approved_path, index=False)
    print(f"Saved: {all_path}")
    print(f"Saved: {approved_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Load step 4-6 outputs
    core_net = cache.load_graph("step4_6_core_net.pkl")
    ranked   = cache.load_df("step4_6_ranked_genes.parquet")

    id_to_symbol     = dict(zip(ranked["string_id"], ranked["symbol"]))
    core_gene_symbols = [id_to_symbol[n] for n in core_net.nodes() if n in id_to_symbol]

    print(f"Querying DGIdb for {len(core_gene_symbols)} genes ...")
    dgi_df      = query_dgidb(core_gene_symbols)
    approved_df = summarise_drug_interactions(dgi_df)

    plot_drug_gene_heatmap(approved_df)
    plot_drugs_per_gene(approved_df)
    save_drug_tables(dgi_df, approved_df)
