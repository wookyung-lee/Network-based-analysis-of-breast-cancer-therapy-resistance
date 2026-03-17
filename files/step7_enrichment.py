"""
Step 7. Enrichment analysis for genes in the core network.

Three databases are queried via Enrichr (gseapy):
  - KEGG_2021_Human           : broad pathway-level view
  - Reactome_2022             : detailed mechanistic pathways
  - GO_Biological_Process_2023: molecular-function vocabulary

Significant terms are defined as adj. p-value < 0.05.
Results are saved to CSV and a bar-chart is generated for KEGG terms.

Loads
-----
.cache/step4_6_core_net.pkl
.cache/step4_6_ranked_genes.parquet   (for id → symbol mapping)

Saves
-----
enrichment_results/*.csv
"""

import os

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cache

DATABASES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "GO_Biological_Process_2023",
]


# ---------------------------------------------------------------------------
# 7.1  Run Enrichr queries
# ---------------------------------------------------------------------------

def run_enrichment(
    gene_symbols: list[str],
    databases: list[str] = DATABASES,
) -> dict[str, pd.DataFrame]:
    """
    Query Enrichr for each database and return a dict of result DataFrames.

    Parameters
    ----------
    gene_symbols : list of HGNC gene symbols
    databases    : list of Enrichr library names

    Returns
    -------
    results : {db_name: DataFrame}
    """
    results = {}
    for db in databases:
        print(f"Querying {db} ...")
        enr = gp.enrichr(
            gene_list=gene_symbols,
            gene_sets=db,
            organism="Human",
            outdir=None,
            verbose=False,
        )
        results[db] = enr.results
        sig = enr.results[enr.results["Adjusted P-value"] < 0.05]
        print(f"  Significant terms (adj. p < 0.05): {len(sig)}")
        print(enr.results[["Term", "Adjusted P-value", "Overlap"]].head(5).to_string(index=False))

    return results


# ---------------------------------------------------------------------------
# 7.2  KEGG bar chart
# ---------------------------------------------------------------------------

def plot_kegg_enrichment(
    kegg_df: pd.DataFrame,
    top_n: int = 20,
    save_path: str | None = "enrichment_KEGG.png",
) -> None:
    """
    Horizontal bar chart of the most significant KEGG pathways.
    Falls back to unadjusted p < 0.1 when no terms pass adj. p < 0.05.
    """
    sig = kegg_df[kegg_df["Adjusted P-value"] < 0.05].head(top_n)

    if sig.empty:
        print("No significant KEGG terms at adj. p < 0.05 — relaxing to p < 0.1")
        sig = kegg_df[kegg_df["P-value"] < 0.1].head(top_n)

    if sig.empty:
        print("No enriched KEGG terms found.")
        return

    sig = sig.copy()
    sig["-log10(Adj.P)"] = -np.log10(sig["Adjusted P-value"].clip(lower=1e-300))
    sig_sorted = sig.sort_values("-log10(Adj.P)")

    fig, ax = plt.subplots(figsize=(9, max(4, len(sig_sorted) * 0.42)))
    ax.barh(sig_sorted["Term"], sig_sorted["-log10(Adj.P)"],
            color="tomato", edgecolor="white")
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--",
               linewidth=0.9, label="$p_{adj} = 0.05$")
    ax.set_xlabel("$-\log_{10}(p_{adj})$")
    # ax.set_title("KEGG Pathway Enrichment — Core Network Genes")
    ax.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# 7.3  Save all results to CSV
# ---------------------------------------------------------------------------

def save_enrichment_results(
    results: dict[str, pd.DataFrame],
    outdir: str = "enrichment_results",
) -> None:
    os.makedirs(outdir, exist_ok=True)
    for db, df in results.items():
        fname = os.path.join(outdir, f"enrichment_{db.replace(' ', '_')}.csv")
        df.to_csv(fname, index=False)
        print(f"Saved: {fname}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Load step 4-6 outputs
    core_net = cache.load_graph("step4_6_core_net.pkl")
    ranked   = cache.load_df("step4_6_ranked_genes.parquet")

    id_to_symbol     = dict(zip(ranked["string_id"], ranked["symbol"]))
    core_gene_symbols = [id_to_symbol[n] for n in core_net.nodes() if n in id_to_symbol]
    print(f"\nGenes submitted to Enrichr: {len(core_gene_symbols)}")

    enrichment_results = run_enrichment(core_gene_symbols)
    plot_kegg_enrichment(enrichment_results["KEGG_2021_Human"])
    save_enrichment_results(enrichment_results)
