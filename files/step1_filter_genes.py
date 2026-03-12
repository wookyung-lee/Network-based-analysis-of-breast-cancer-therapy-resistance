"""
Step 1. Prune the gene list by removing all genes that do not have
expression above 1 TPM in at least 3 samples.

TPM (Transcripts Per Million): a normalized measure of gene expression in
RNA-Seq data. It estimates how abundant a gene's transcripts are relative
to all transcripts in the sample.

Saves
-----
.cache/step1_filtered_genes.parquet
"""

import pandas as pd

import cache


def filter_by_tpm(filepath: str, min_tpm: float = 1.0, min_samples: int = 3) -> pd.DataFrame:
    """
    Load a TPM matrix and keep only genes expressed above `min_tpm`
    in at least `min_samples` samples.

    Parameters
    ----------
    filepath    : path to the TSV file (first column must be GeneID)
    min_tpm     : expression threshold (default 1 TPM)
    min_samples : minimum number of samples that must exceed the threshold

    Returns
    -------
    filtered_genes : DataFrame with passing genes
    """
    tpm = pd.read_csv(filepath, sep="\t")

    # iloc[:, 1:] skips the GeneID column when evaluating the threshold
    mask = (tpm.iloc[:, 1:] > min_tpm).sum(axis=1) >= min_samples
    filtered_genes = tpm.loc[mask].copy()

    print(f"Genes after TPM filter: {filtered_genes.shape[0]}")
    return filtered_genes


if __name__ == "__main__":
    filtered_genes = filter_by_tpm(
        "other_dataset/GSE162187_raw_counts_GRCh38.p13_NCBI.tsv"
    )
    cache.save_df(filtered_genes, "step1_filtered_genes.parquet")
    print(filtered_genes.head())
