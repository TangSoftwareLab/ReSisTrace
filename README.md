# ReSisTrace
The data analysis codes for the figures of ReSisTrace manuscript.

# Input data

"./data/input/all_bt_sample_merge_filtered_gene_sisters.RDS" is a Seurat object. It merges the pre-treatment samples from all of the 8 experiments (4 treatments, 2 replicates). Lowly expressed genes whose total normalised expressions are less than 1.25 (10% quantile) were removed from this object and the sister cells with the Euclidean distance larger than 14.9 (90% quantile) were not marked as sisters.