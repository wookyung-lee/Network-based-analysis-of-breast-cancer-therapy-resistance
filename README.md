# Network-based-analysis-of-melanoma-therapy-resistance
Project from the Course "Introduction to simulation, network and data analysis in Medical Systems Biology (WS2526)"

Datasets: \\

1. Dataset relevant for step 1 (TPM):
GSE78220_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
link: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE78220 
structure of tsv file:
- Rows = genes
- First column = GeneID
- Remaining columns = samples (GSM IDs)
- Values = TPM

2. Dataset relevant for step 2 (protein links):
9606.protein.links.v12.0.txt
link: https://string-db.org/cgi/download?sessionId=bOLP4AG5dmGT&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt

3. Dataset relevant for step 2 (protein_id - symbol (preferred name)):
9606.protein.info.v12.0.txt
This is used to map STRING protein IDs (used in Dataset 2) to Gene symbols (used in Dataset 1)
link: https://string-db.org/cgi/download?sessionId=bOLP4AG5dmGT&species_text=Homo+sapiens&settings_expanded=0&min_download_score=0&filter_redundant_pairs=0&delimiter_type=txt

4. Dataset relevant for step 3 (raw counts per sample): 
GSE78220_raw_counts_GRCh38.p13_NCBI.tsv
link: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE78220

5. Dataset relevant for step 3 (PD, CR response):
GSE78220_series_matrix.txt
link: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE78220

6. 

These datasets are from paper: \\
Hugo et al. 2016. PubMed identifier PMID: 26997480
link: https://pubmed.ncbi.nlm.nih.gov/26997480/ 