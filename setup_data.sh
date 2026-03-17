#!/bin/bash
set -euo pipefail
mkdir -p other_dataset step2
echo "--- Starting Data Download ---"
# -- 1. GSE162187 TPM counts ---------------------------------------------------
echo ""
echo "Fetching GSE162187 TPM counts..."
curl -L -o other_dataset/GSE162187_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz \
    "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE162187&format=file&file=GSE162187_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
# -- 2. GSE162187 raw counts ---------------------------------------------------
echo ""
echo "Fetching GSE162187 raw counts..."
curl -L -o other_dataset/GSE162187_raw_counts_GRCh38.p13_NCBI.tsv.gz \
    "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE162187&format=file&file=GSE162187_raw_counts_GRCh38.p13_NCBI.tsv.gz"
# -- 3. GSE162187 series matrix (sample metadata) ------------------------------
echo ""
echo "Fetching GSE162187 series matrix..."
curl -L -o other_dataset/GSE162187_series_matrix.txt.gz \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162187/matrix/GSE162187_series_matrix.txt.gz"
# -- 5. Decompress GEO files ---------------------------------------------------
echo ""
echo "Decompressing GEO files..."
for f in other_dataset/*.gz; do
    [ -f "$f" ] && gunzip -v "$f"
done
# -- 6. STRING database v12.0 (human, taxon 9606) ------------------------------
echo ""
echo "Fetching STRING v12.0 files..."
curl -L -o other_dataset/9606.protein.links.v12.0.txt.gz \
    "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
curl -L -o other_dataset/9606.protein.info.v12.0.txt.gz \
    "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
echo "Decompressing STRING files..."
gunzip other_dataset/9606.protein.links.v12.0.txt.gz
gunzip other_dataset/9606.protein.info.v12.0.txt.gz
# -- Summary -------------------------------------------------------------------
echo ""
echo "--- Setup Complete ---"
echo ""
echo "other_dataset/ contains:"
ls -lh other_dataset/