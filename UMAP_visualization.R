library(Seurat)
library(irlba)
library(umap)
library(plyr)

# Transformation ---------------------------------------------------------------

data = readRDS("./data/input/merge_all_sample_with_CNV_subclone.RDS")
data = FindVariableFeatures(data)
data = ScaleData(data)
saveRDS(data, file = "./data/input/merged_KURAMOCHI_samples1_scaled_new.RDS")
print("transformation ok!")

# Visualization ----------------------------------------------------------------

data = readRDS("./data/input/merged_KURAMOCHI_samples1_scaled_new.RDS")
data2 = data@assays$RNA@scale.data
dim(data2)

print("run pca")
P <- prcomp_irlba(data2, n = 50, scale. = T)
summary(P)

print("run umap")
data2.umap = umap(P$rotation)
saveRDS(data2.umap, file = "./data/input/merged_KURAMOCHI_samples1_SCT_umap_2000_n50_scaleT.RDS")

print("visualization ok!")

