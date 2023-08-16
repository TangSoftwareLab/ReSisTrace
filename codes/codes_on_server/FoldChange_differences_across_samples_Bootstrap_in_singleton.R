library(openxlsx)
library(tidyverse)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)

b <- as.numeric(args[1]) # number of bootstrap samples
cat("Number of bootstrap samples: ", b, "\n")

output_dir <- paste0("../data/output/FoldChange_differences_across_samples/bootstrap_in_singleton/", b, "_iterations/")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# load data
## scRNA-seq data (treatment data)
scRNA <- readRDS("../data/raw/all_bt_sample_merge_filtered_gene_sisters.RDS")
meta <- scRNA@misc$gene_meta
## scRNA-seq data (growth control)
control <- readRDS("../data/output/scLT_ctrl_1/preprocess/bt_sample_ready.RDS")
control$sisters_sample <- paste0("ctrl_1_", control$sisters)
control$sample <- "ctrl_1"
control$true_sister <- NA
control <- RenameCells(control, new.names = paste0("ctrl_1__", colnames(control)))
control <- subset(control, features = rownames(scRNA))
scRNA <- merge(scRNA, control)

control <- readRDS("../data/output/scLT_ctrl_2/preprocess/bt_sample_ready.RDS")
control$sisters_sample <- paste0("ctrl_2_", control$sisters)
control$sample <- "ctrl_2"
control$true_sister <- NA
control <- RenameCells(control, new.names = paste0("ctrl_2__", colnames(control)))
control <- subset(control, features = rownames(scRNA))
scRNA <- merge(scRNA, control)
scRNA@misc$gene_meta <- meta
rm(control)

# bootstrapping
# https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r

samples <- c("Olaparib", "Carbo", "NK", "ctrl")
fc_table <- NULL

for (s in samples){
  pres <- subset(
    scRNA, 
    cells = colnames(scRNA)[scRNA@meta.data$drugSens == "pre-sensitive" & 
                              !is.na(scRNA@meta.data$drugSens) &
			      is.na(scRNA@meta.data$sisters) &
                              scRNA@meta.data$sample %in% paste0(s, c("_1", "_2"))]
  )
  pres <- pres@assays$RNA@data
  prer <- subset(
    scRNA, 
    cells = colnames(scRNA)[scRNA@meta.data$drugSens == "pre-resistant" & 
                              !is.na(scRNA@meta.data$drugSens) &
			      is.na(scRNA@meta.data$sisters) &
                              scRNA@meta.data$sample %in% paste0(s, c("_1", "_2"))]
  )
  prer <- prer@assays$RNA@data
  
  fc <- NULL
  for(i in 1:b){
    pres_mean_exp <- log2(
      rowMeans(
        expm1(
          pres[, sample(colnames(pres), ncol(pres), replace = T)]
        ),
        na.rm = TRUE
      ) + 1
    )
    prer_mean_exp <- log2(
      rowMeans(
        expm1(
          prer[, sample(colnames(prer), ncol(prer), replace = T)]
        ),
        na.rm = TRUE
      ) + 1
    )
    fc <- cbind(fc, prer_mean_exp - pres_mean_exp)
  }
  fc <- data.frame(fc)
  colnames(fc) <- paste0(s, 1:b)
  fc$gene <- rownames(fc)
  if (is.null(fc_table)){
    fc_table <- fc
  } else {
    fc_table <- full_join(fc_table, fc, by = "gene")
  }
}
fc_table <- select(fc_table, gene, everything())

fc_table <- left_join(
  fc_table,
  scRNA@misc$gene_meta,
  by = c("gene" = "ensembl")
)
write.csv(fc_table, paste0(output_dir, "/bootstrap_fc_table.csv"), row.names = FALSE)

combo <- combn(samples, 2)

stat_table <- NULL
for (i in 1:ncol(combo)){
  # confidence interval of fc1-fc2
  tmp <- fc_table[, startsWith(colnames(fc_table), combo[1, i])] -
    fc_table[, startsWith(colnames(fc_table), combo[2, i])]
  ci <- apply(tmp, 1, function(x){quantile(x,c(0.025,0.975))})
  ci <- t(ci)
  # p-value
  z <- abs(rowMeans(tmp, na.rm = TRUE)/sqrt(rowSums((tmp-rowMeans(tmp, na.rm = TRUE))^2)/(b-1)))
  p_value <- exp(-0.717*z-0.416*z*z)
  tmp_stat <- cbind(ci, z, p_value)
  colnames(tmp_stat) <- paste(paste0(combo[, i], collapse = "_"), colnames(tmp_stat))
  tmp_stat <- as.data.frame(tmp_stat)
  tmp_stat$gene <- fc_table$gene
  if (is.null(stat_table)){
    stat_table <- tmp_stat
  } else {
    stat_table <- full_join(stat_table, tmp_stat, by = "gene")
  }
}

stat_table <- left_join(
  stat_table,
  scRNA@misc$gene_meta,
  by = c("gene" = "ensembl")
)
write.csv(stat_table, paste0(output_dir, "bootstrap_stat_table.csv"), row.names = FALSE)
