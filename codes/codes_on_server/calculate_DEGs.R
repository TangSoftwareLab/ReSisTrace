library(Seurat)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
s <- args[1]
options(future.globals.maxSize= 1024*1024^300)

data <- readRDS(paste0("../data/output/", s, "/preprocess/bt_sample_ready.RDS"))
data
Idents(data) <- "drugSens"
output_dir <- paste0("../data/output/", s, "/DEG/")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

low_exp <- read.csv("../data/raw/KURAMOCHI_filtered_low_exp_genes.csv")
data <- subset(data, features = setdiff(rownames(data), low_exp$gene))

# preR-preS in singleton
print("preR-preS in singleton")
deg_prer_pres <- FindMarkers(
    subset(data, cells = colnames(data)[is.na(data$sisters)]),
    ident.1 = "pre-resistant", 
    ident.2 = "pre-sensitive",
    test.use = "wilcox", # old version is "wilcox"
    logfc.threshold = 0,
    min.pct = 0,
    verbose = TRUE
  )
deg_prer_pres <- deg_prer_pres %>% 
  rownames_to_column(var = "ensembl") %>% 
  left_join(as.data.frame(data@misc$gene_meta), by = "ensembl")

write.csv(deg_prer_pres, paste0(output_dir, "/deg_prer_pres_in_singleton_wilcox.csv"), row.names = FALSE)