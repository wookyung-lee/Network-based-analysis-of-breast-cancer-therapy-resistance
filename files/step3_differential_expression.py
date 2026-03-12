"""
Step 3. Compute p-values and log2 fold-change for every surviving gene.
Apply Benjamini-Hochberg (BH) correction, then generate a volcano plot
and a PCA plot.

Definitions
-----------
Fold Change (FC)   : mean(group2) / mean(group1) — how much expression differs.
log2FC             : log2(FC); positive = higher in Sensitive, negative = higher in Resistant.
BH correction      : adjusts p-values for the number of simultaneous tests to
                     control the False Discovery Rate (FDR).
Volcano plot       : scatter of log2FC (x) vs −log10(p-value) (y).
PCA                : principal component analysis — checks whether samples
                     cluster by response group.

Loads
-----
.cache/step2_filtered_genes.parquet

Saves
-----
.cache/step3_filtered_genes.parquet   (with log2FC, pval, pval_adj columns added)
.cache/step3_counts.parquet           (expression matrix used for testing)
.cache/step3_sample_groups.json       (resistant_gsm and sensitive_gsm lists)
"""

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

import cache

EPS = 1e-6


# ---------------------------------------------------------------------------
# 3.1  Parse sample metadata from a GEO series matrix
# ---------------------------------------------------------------------------

def parse_series_matrix(filepath: str) -> tuple[dict, list, list]:
    """
    Extract GSM IDs and treatment-response labels from a GEO series matrix.

    Returns
    -------
    gsm_group_info  : {gsm_id: response_label}
    resistant_gsm   : list of GSM IDs labelled "Resistant"
    sensitive_gsm   : list of GSM IDs labelled "Sensitive"
    """
    with open(filepath) as f:
        lines = f.readlines()

    char_lines = [l.strip() for l in lines if l.startswith("!Sample_characteristics_ch1")]
    title_line = [l.strip() for l in lines if l.startswith("!Sample_title")][0]
    gsm_line   = [l.strip() for l in lines if l.startswith("!Sample_geo_accession")][0]

    response_row = [l for l in char_lines if "treatment response:" in l][0]

    sample_titles  = title_line.split("\t")[1:]
    response_data  = response_row.split("\t")[1:]
    gsm_ids        = [x.replace('"', '') for x in gsm_line.split("\t")[1:]]

    patients  = [x.replace('"', '').split('_')[-1] for x in sample_titles]
    responses = [
        re.search(r"treatment response: (.+?)$", x.replace('"', '')).group(1)
        for x in response_data
    ]

    gsm_group_info = dict(zip(gsm_ids, responses))
    resistant_gsm  = [g for g, r in gsm_group_info.items() if r == "Resistant"]
    sensitive_gsm  = [g for g, r in gsm_group_info.items() if r == "Sensitive"]

    print(f"Resistant: {len(resistant_gsm)}  |  Sensitive: {len(sensitive_gsm)}")
    return gsm_group_info, resistant_gsm, sensitive_gsm


# ---------------------------------------------------------------------------
# 3.2  Differential expression — p-values and log2FC
# ---------------------------------------------------------------------------

def compute_differential_expression(
    counts_path: str,
    filtered_genes: pd.DataFrame,
    resistant_gsm: list,
    sensitive_gsm: list,
) -> pd.DataFrame:
    """
    For every gene that survived the TPM filter, run a Mann-Whitney U test
    (Sensitive vs Resistant) and calculate log2FC.  Applies BH correction.

    Returns
    -------
    filtered_genes : input DataFrame extended with log2FC, pval, pval_adj
    counts         : expression matrix (genes × samples) used for the test
    """
    counts = pd.read_csv(counts_path, sep="\t")
    counts = counts[counts["GeneID"].isin(filtered_genes["GeneID"])].set_index("GeneID")

    res_cols = [c for c in counts.columns if c in resistant_gsm]
    sen_cols = [c for c in counts.columns if c in sensitive_gsm]
    print(f"Resistant columns: {len(res_cols)}  |  Sensitive columns: {len(sen_cols)}")

    pvals, log2fc = [], []

    for gene in counts.index:
        res_expr = counts.loc[gene, res_cols].values.astype(float)
        sen_expr = counts.loc[gene, sen_cols].values.astype(float)

        fc = np.log2((np.mean(sen_expr) + EPS) / (np.mean(res_expr) + EPS))
        log2fc.append(fc)

        try:
            _, p = mannwhitneyu(sen_expr, res_expr, alternative="two-sided")
        except ValueError:
            p = 1.0
        pvals.append(p)

    _, pvals_adj, _, _ = multipletests(pvals, method="fdr_bh")

    stat_df = pd.DataFrame({
        "GeneID"  : counts.index.astype(str),
        "log2FC"  : log2fc,
        "pval"    : pvals,
        "pval_adj": pvals_adj,
    }).reset_index(drop=True)

    filtered_genes = filtered_genes.copy()
    filtered_genes["GeneID"] = filtered_genes["GeneID"].astype(str)
    filtered_genes = filtered_genes.merge(stat_df, on="GeneID", how="left")

    print(f"Genes tested: {len(stat_df)}")
    print(f"Raw p < 0.05: {(stat_df['pval'] < 0.05).sum()}")
    print(f"Adj p < 0.05: {(stat_df['pval_adj'] < 0.05).sum()}")

    return filtered_genes, counts


# ---------------------------------------------------------------------------
# 3.3  Volcano plot
# ---------------------------------------------------------------------------

def plot_volcano(
    filtered_genes: pd.DataFrame,
    use_adjusted: bool = False,
    fc_thresh: float = 1.0,
    p_thresh: float = 0.05,
    save_path: str | None = None,
) -> None:
    """
    Draw a volcano plot.

    Parameters
    ----------
    use_adjusted : if True, plot adjusted p-values on the y-axis
    fc_thresh    : |log2FC| threshold for colouring (default 1)
    p_thresh     : p-value threshold for colouring (default 0.05)
    save_path    : file path to save the figure; None = do not save
    """
    p_col  = "pval_adj" if use_adjusted else "pval"
    p_label = "adj. p-value" if use_adjusted else "p-value"

    y = -np.log10(filtered_genes[p_col] + EPS)

    colors = []
    for fc, p in zip(filtered_genes["log2FC"], filtered_genes[p_col]):
        if p < p_thresh and fc >  fc_thresh:
            colors.append("red")     # up in Sensitive
        elif p < p_thresh and fc < -fc_thresh:
            colors.append("blue")    # up in Resistant
        else:
            colors.append("grey")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(filtered_genes["log2FC"], y, c=colors, alpha=0.5, s=10)
    ax.axhline(-np.log10(p_thresh),  color="black", linestyle="--", linewidth=0.8)
    ax.axvline( fc_thresh,           color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_thresh,           color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change (Sensitive vs Resistant)")
    ax.set_ylabel(f"-log10({p_label})")
    ax.set_title(f"Volcano Plot: Sensitive vs Resistant"
                 + (" (FDR-adjusted)" if use_adjusted else ""))
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# 3.4  PCA
# ---------------------------------------------------------------------------

def plot_pca(
    counts: pd.DataFrame,
    group_cols: dict,          # e.g. {"Resistant": [...], "Sensitive": [...]}
    group_colors: dict | None = None,
    save_path: str | None = None,
) -> None:
    """
    Run PCA on the supplied expression matrix and plot PC1 vs PC2.

    Parameters
    ----------
    counts       : gene × sample expression DataFrame (index = GeneID)
    group_cols   : {label: [sample_col_names]}
    group_colors : {label: colour_string}; auto-assigned if None
    save_path    : optional output path
    """
    default_colors = ["steelblue", "tomato", "forestgreen", "orange", "purple"]
    if group_colors is None:
        group_colors = {
            label: default_colors[i % len(default_colors)]
            for i, label in enumerate(group_cols)
        }

    all_cols, labels = [], []
    for label, cols in group_cols.items():
        valid = [c for c in cols if c in counts.columns]
        all_cols.extend(valid)
        labels.extend([label] * len(valid))

    pca_data = counts[all_cols].T
    pca_log  = np.log2(pca_data + 1)
    scaled   = StandardScaler().fit_transform(pca_log)

    pca        = PCA(n_components=2)
    components = pca.fit_transform(scaled)

    fig, ax = plt.subplots(figsize=(7, 5))
    for label in group_cols:
        idx = [i for i, l in enumerate(labels) if l == label]
        ax.scatter(
            components[idx, 0], components[idx, 1],
            label=label, color=group_colors[label], s=60,
        )
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA: " + " vs ".join(group_cols.keys()))
    ax.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Load step 2 output
    filtered_genes = cache.load_df("step2_filtered_genes.parquet")

    _, resistant_gsm, sensitive_gsm = parse_series_matrix(
        "other_dataset/GSE162187_series_matrix.txt"
    )

    filtered_genes, counts = compute_differential_expression(
        counts_path="other_dataset/GSE162187_raw_counts_GRCh38.p13_NCBI.tsv",
        filtered_genes=filtered_genes,
        resistant_gsm=resistant_gsm,
        sensitive_gsm=sensitive_gsm,
    )

    # Save step 3 outputs
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
