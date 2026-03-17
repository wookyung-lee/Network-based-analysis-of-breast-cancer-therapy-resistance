"""
main.py — Run the full pipeline end-to-end.

Each step saves its outputs to .cache/ so you can also run the steps
individually:

    python step1_filter_genes.py
    python step2a_network_pruning.py 
    python step2b_network_metrics.py
    python step3_differential_expression.py
    python step4_6_scoring_network.py
    python step7_enrichment.py
    python step8_drug_interactions.py

Edit the paths in CONFIG before running.
"""

import cache
from step1_filter_genes import filter_by_tpm
from step2a_network_pruning import build_pruned_network, map_entrez_to_symbols
from step2b_network_metrics import compute_network_metrics
from step3_differential_expression import (
    compute_differential_expression,
    parse_series_matrix,
    plot_pca,
    plot_volcano,
)
from step4_6_scoring_network import (
    build_core_network,
    plot_core_network,
    rank_genes,
    score_genes,
)
from step7_enrichment import plot_kegg_enrichment, run_enrichment, save_enrichment_results
from step8_drug_interactions import (
    plot_drug_gene_heatmap,
    plot_drugs_per_gene,
    query_dgidb,
    save_drug_tables,
    summarise_drug_interactions,
)

# ---------------------------------------------------------------------------
# Configuration — edit these paths to match your local file layout
# ---------------------------------------------------------------------------
CONFIG = {
    "raw_counts"   : "other_dataset/GSE162187_norm_counts_TPM_GRCh38.p13_NCBI.tsv",
    "series_matrix": "other_dataset/GSE162187_series_matrix.txt",
    "string_links" : "other_dataset/9606.protein.links.v12.0.txt",
    "string_info"  : "other_dataset/9606.protein.info.v12.0.txt",
}


def main():
    # ------------------------------------------------------------------
    # Step 1
    # ------------------------------------------------------------------
    print("=" * 60)
    print("STEP 1 — TPM filtering")
    print("=" * 60)
    filtered_genes = filter_by_tpm(CONFIG["raw_counts"])
    cache.save_df(filtered_genes, "step1_filtered_genes.parquet")

    # ------------------------------------------------------------------
    # Step 2a
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 2a — Network construction")
    print("=" * 60)
    filtered_genes = map_entrez_to_symbols(filtered_genes)
    filtered_genes, G = build_pruned_network(
        filtered_genes,
        links_path=CONFIG["string_links"],
        info_path=CONFIG["string_info"],
    )
    cache.save_df(filtered_genes, "step2a_filtered_genes.parquet")
    cache.save_graph(G, "step2a_graph.pkl")

    # ------------------------------------------------------------------
    # Step 2b
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 2b — ND / BC")
    print("=" * 60)
    filtered_genes = compute_network_metrics(filtered_genes, G)
    cache.save_df(filtered_genes, "step2b_filtered_genes.parquet")

    # ------------------------------------------------------------------
    # Step 3
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 3 — Differential expression, volcano plot, PCA")
    print("=" * 60)
    _, resistant_gsm, sensitive_gsm = parse_series_matrix(CONFIG["series_matrix"])

    filtered_genes, counts = compute_differential_expression(
        counts_path=CONFIG["raw_counts"],
        filtered_genes=filtered_genes,
        resistant_gsm=resistant_gsm,
        sensitive_gsm=sensitive_gsm,
    )
    cache.save_df(filtered_genes, "step3_filtered_genes.parquet")
    cache.save_df(counts, "step3_counts.parquet")
    cache.save_json(
        {"resistant": resistant_gsm, "sensitive": sensitive_gsm},
        "step3_sample_groups.json",
    )

    plot_volcano(filtered_genes, use_adjusted=False, save_path="volcano_raw.png")
    plot_volcano(filtered_genes, use_adjusted=True,  save_path="volcano_adj.png")
    plot_pca(
        counts,
        group_cols={"Resistant": resistant_gsm, "Sensitive": sensitive_gsm},
        save_path="pca_plot.png",
    )

    # ------------------------------------------------------------------
    # Steps 4 – 6
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEPS 4–6 — Scoring, ranking, core network")
    print("=" * 60)
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
    plot_core_network(core_net, id_to_symbol, id_to_log2fc, id_to_score)

    # ------------------------------------------------------------------
    # Step 7
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 7 — Enrichment analysis")
    print("=" * 60)
    core_gene_symbols = [id_to_symbol[n] for n in core_net.nodes() if n in id_to_symbol]
    print(f"Genes submitted to Enrichr: {len(core_gene_symbols)}")

    enrichment_results = run_enrichment(core_gene_symbols)
    plot_kegg_enrichment(enrichment_results["KEGG_2021_Human"])
    save_enrichment_results(enrichment_results)

    # ------------------------------------------------------------------
    # Step 8
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 8 — Drug-gene interactions (DGIdb)")
    print("=" * 60)
    dgi_df      = query_dgidb(core_gene_symbols)
    approved_df = summarise_drug_interactions(dgi_df)

    plot_drug_gene_heatmap(approved_df)
    plot_drugs_per_gene(approved_df)
    save_drug_tables(dgi_df, approved_df)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
