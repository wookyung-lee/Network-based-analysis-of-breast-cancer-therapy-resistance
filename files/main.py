import os
import cache

# 1. Imports from your files (only what is explicitly defined as a function)
from step1_filter_genes import filter_by_tpm
from step2a_network_pruning import map_entrez_to_symbols, build_pruned_network
from step2b_network_metrics import compute_network_metrics
from step3_differential_expression import (
    parse_series_matrix, compute_differential_expression, plot_volcano, plot_pca
)
from step4_6_scoring_network import (
    score_genes, rank_genes, build_core_network, plot_core_network
)
from step7_enrichment import run_enrichment, save_enrichment_results
from step8_drug_interactions import (
    query_dgidb, summarise_drug_interactions, save_drug_tables
)

def run_pipeline():
    # --- Configuration ---
    RAW_COUNTS = "other_dataset/GSE162187_raw_counts_GRCh38.p13_NCBI.tsv"
    SERIES_MAT = "other_dataset/GSE162187_series_matrix.txt"
    STRING_LINKS = "other_dataset/9606.protein.links.v12.0.txt"
    STRING_INFO = "other_dataset/9606.protein.info.v12.0.txt"
    
    OUT_DIR = "pipeline_outputs"
    os.makedirs(OUT_DIR, exist_ok=True)

    print("🚀 Starting End-to-End Pipeline")

    # --- Step 1 & 2: Filtering and Network ---
    df = filter_by_tpm(RAW_COUNTS)
    df = map_entrez_to_symbols(df)
    df, G = build_pruned_network(df, STRING_LINKS, STRING_INFO)
    df = compute_network_metrics(df, G)
    cache.save_df(df, "step2_final.parquet")
    cache.save_graph(G, "step2_graph.pkl")

    # --- Step 3: DE Analysis & Figures ---
    _, res_gsm, sen_gsm = parse_series_matrix(SERIES_MAT)
    df, counts = compute_differential_expression(RAW_COUNTS, df, res_gsm, sen_gsm)
    
    # These functions ARE defined in your step3 file
    plot_volcano(df, save_path=f"{OUT_DIR}/01_volcano.png")
    
    plot_pca(counts, {"Resistant": res_gsm, "Sensitive": sen_gsm}, save_path=f"{OUT_DIR}/02_pca.png")
    cache.save_df(df, "step3_de_results.parquet")

    # --- Step 4-6: Scoring & Core Network ---
    df = score_genes(df)
    ranked = rank_genes(df)
    core_net, _, _ = build_core_network(ranked, G)
    
    id_to_symbol = dict(zip(ranked["string_id"], ranked["symbol"]))
    id_to_score = dict(zip(ranked["string_id"], ranked["score"]))
    
    # This function IS defined in your step4_6 file
    plot_core_network(core_net, id_to_symbol, id_to_score, save_path=f"{OUT_DIR}/03_network.png")
    
    
    cache.save_df(ranked, "step4_6_ranked.parquet")
    cache.save_graph(core_net, "step4_6_core.pkl")

    # --- Step 7: Enrichment ---
    core_symbols = [id_to_symbol[n] for n in core_net.nodes() if n in id_to_symbol]
    enr_results = run_enrichment(core_symbols)
    save_enrichment_results(enr_results, outdir=f"{OUT_DIR}/enrichment")
    
    # Note: I removed plot_kegg because it is hidden inside your 'if __name__ == "__main__"' block.
    print("[Notice] Enrichment plots skipped because function is not exported in step7.")
    

    # --- Step 8: Drugs ---
    dgi_df = query_dgidb(core_symbols)
    app_df = summarise_drug_interactions(dgi_df)
    save_drug_tables(dgi_df, app_df, f"{OUT_DIR}/05_all_drugs.csv", f"{OUT_DIR}/06_app_drugs.csv")
    
    # This function IS defined in your step8 file
    # plot_drug_interactions(app_df, save_path=f"{OUT_DIR}/07_drugs.png")

    print(f"\n✅ Pipeline Complete. All available figures and data saved to {OUT_DIR}")

if __name__ == "__main__":
    run_pipeline()
