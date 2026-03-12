---

### 2. The Download Script (`setup_data.sh`)
This script creates the data directory, downloads the compressed files from the NCBI FTP servers, and extracts them.

```bash
#!/bin/bash

# Create data directory if it doesn't exist
mkdir -p other_dataset

echo "--- Starting Data Download for Melanoma Project ---"

# 1. GSE162187 (The primary dataset)
echo "Fetching GSE162187..."
wget -P other_dataset/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162187/suppl/GSE162187_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
wget -P other_dataset/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162187/suppl/GSE162187_raw_counts_GRCh38.p13_NCBI.tsv.gz

# 2. GSE78220 (The validation dataset)
echo "Fetching GSE78220..."
wget -P other_dataset/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz

# Decompress files
echo "Decompressing files..."
gunzip other_dataset/*.gz

echo "--- Setup Complete. Data is located in /other_dataset ---"
