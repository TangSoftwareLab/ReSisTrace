library(tidyverse)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
print(args)
s <- args[1] # CNV subclone
print(s)
output_dir <- "../data/output/CNV_DEG/"
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
data <- readRDS("../data/raw/Jing_code_data/merged_KURAMOCHI_samples1_scaled_new.RDS")

subclone <- unique(data$CNV_subclone)
sub_data <- subset(data, cells = colnames(data)[data$time_point == "BeforeTreatment"])

Idents(sub_data) <- "CNV_subclone"

deg <- FindMarkers(
  sub_data,
  ident.1 = s,
  ident.2 = NULL,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  return.thresh = 1,
  only.pos = FALSE
)
write.csv(deg, paste0(output_dir, s, "_vs_others_DEG.csv"))
