"""
Step 7. Enrichment analysis — two gene sets:

  A) Core network genes (from Steps 4-6 scoring pipeline)
  B) Significant DEGs from the volcano plot (padj < 0.05, |log2FC| > 1),
     split further into up-in-Sensitive and up-in-Resistant.

Three databases are queried via Enrichr (gseapy) for each gene set:
  - KEGG_2021_Human
  - Reactome_2022
  - GO_Biological_Process_2023

Loads
-----
.cache/step4_6_core_net.pkl
.cache/step4_6_ranked_genes.parquet
.cache/step3_sig_genes.parquet

Saves
-----
enrichment_results/                  <- core network
enrichment_results_up_sensitive/     <- volcano: up in Sensitive
enrichment_results_up_resistant/     <- volcano: up in Resistant
enrichment_results_all_sig/          <- volcano: all significant
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
    """
    results = {}
    for db in databases:
        print(f"  Querying {db} ...")
        enr = gp.enrichr(
            gene_list=gene_symbols,
            gene_sets=db,
            organism="Human",
            outdir=None,
            verbose=False,
        )
        results[db] = enr.results
        sig = enr.results[enr.results["Adjusted P-value"] < 0.05]
        print(f"    Significant terms (adj. p < 0.05): {len(sig)}")
        print(enr.results[["Term", "Adjusted P-value", "Overlap"]].head(5).to_string(index=False))

    return results


# ---------------------------------------------------------------------------
# 7.2  KEGG bar chart
# ---------------------------------------------------------------------------

def plot_kegg_enrichment(
    kegg_df: pd.DataFrame,
    title: str = "KEGG Pathway Enrichment",
    top_n: int = 20,
    save_path: str | None = "enrichment_KEGG.png",
) -> None:
    """
    Horizontal bar chart of the most significant KEGG pathways.
    Falls back to unadjusted p < 0.1 when no terms pass adj. p < 0.05.
    """
    sig = kegg_df[kegg_df["Adjusted P-value"] < 0.05].head(top_n)

    if sig.empty:
        print(f"  No significant KEGG terms at adj. p < 0.05 — relaxing to p < 0.1")
        sig = kegg_df[kegg_df["P-value"] < 0.1].head(top_n)

    if sig.empty:
        print(f"  No enriched KEGG terms found for: {title}")
        return

    sig = sig.copy()
    sig["-log10(Adj.P)"] = -np.log10(sig["Adjusted P-value"].clip(lower=1e-300))
    sig_sorted = sig.sort_values("-log10(Adj.P)")

    fig, ax = plt.subplots(figsize=(9, max(4, len(sig_sorted) * 0.42)))
    ax.barh(sig_sorted["Term"], sig_sorted["-log10(Adj.P)"],
            color="tomato", edgecolor="white")
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--",
               linewidth=0.9, label="adj. p = 0.05")
    ax.set_xlabel("-log10(Adjusted P-value)")
    # ax.set_title(title)
    ax.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
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
        print(f"  Saved: {fname}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    # ── A. Core network enrichment ────────────────────────────────────────────
    print("=" * 60)
    print("A. Enrichment: Core Network Genes")
    print("=" * 60)

    core_net = cache.load_graph("step4_6_core_net.pkl")
    ranked   = cache.load_df("step4_6_ranked_genes.parquet")

    id_to_symbol      = dict(zip(ranked["string_id"], ranked["symbol"]))
    core_gene_symbols = [id_to_symbol[n] for n in core_net.nodes() if n in id_to_symbol]
    print(f"Genes submitted: {len(core_gene_symbols)}")

    core_results = run_enrichment(core_gene_symbols)
    plot_kegg_enrichment(
        core_results["KEGG_2021_Human"],
        title="KEGG Enrichment — Core Network Genes",
        save_path="enrichment_KEGG_core.png",
    )
    save_enrichment_results(core_results, outdir="enrichment_results_core")

    # ── B. Volcano significant genes enrichment ───────────────────────────────
    print("\n" + "=" * 60)
    print("B. Enrichment: Significant Volcano DEGs")
    print("=" * 60)

    sig_genes    = cache.load_df("step3_sig_genes.parquet")
    up_sensitive = sig_genes[sig_genes["direction"] == "up_in_sensitive"]["symbol"].dropna().tolist()
    up_resistant = sig_genes[sig_genes["direction"] == "up_in_resistant"]["symbol"].dropna().tolist()
    all_sig      = sig_genes["symbol"].dropna().tolist()

    print(f"Up in Sensitive : {len(up_sensitive)} genes")
    print(f"Up in Resistant : {len(up_resistant)} genes")
    print(f"All significant : {len(all_sig)} genes")

    gene_sets = {
        "up_sensitive": up_sensitive,
        "up_resistant": up_resistant,
        "all_sig"     : all_sig,
    }
    titles = {
        "up_sensitive": "KEGG Enrichment — Up in Sensitive (red)",
        "up_resistant": "KEGG Enrichment — Up in Resistant (blue)",
        "all_sig"     : "KEGG Enrichment — All Significant DEGs",
    }

    for label, genes in gene_sets.items():
        if not genes:
            print(f"\nSkipping {label} — no genes.")
            continue

        print(f"\n--- {label} ({len(genes)} genes) ---")
        results = run_enrichment(genes)
        plot_kegg_enrichment(
            results["KEGG_2021_Human"],
            title=titles[label],
            save_path=f"enrichment_KEGG_{label}.png",
        )
        save_enrichment_results(results, outdir=f"enrichment_results_{label}")
