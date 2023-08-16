library(openxlsx)
library(tidyverse)
library(Seurat)
library(purrr)
library(modelr)
library(parallel)

args <- commandArgs(trailingOnly=TRUE)

p <- as.numeric(args[1]) # number of permutation samples
cat("Number of permutation samples: ", p, "\n")

output_dir <- paste0("../data/output/FoldChange_differences_across_samples/permutation_in_singleton/", p, "_iterations/")
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
scRNA <- merge(scRNA, control)

control <- readRDS("../data/output/scLT_ctrl_2/preprocess/bt_sample_ready.RDS")
control$sisters_sample <- paste0("ctrl_2_", control$sisters)
control$sample <- "ctrl_2"
control$true_sister <- NA
control <- RenameCells(control, new.names = paste0("ctrl_2__", colnames(control)))
scRNA <- merge(scRNA, control)
scRNA@misc$gene_meta <- meta
rm(control)
scRNA <- subset(scRNA, cells = colnames(scRNA)[is.na(scRNA@meta.data$sisters)])

# Extract classes

prer <- scRNA@meta.data %>% 
  rownames_to_column("cell") %>% 
  select(cell, drugSens, sample) %>% 
  filter(
    drugSens == "pre-resistant" &
      !is.na(drugSens)
  )
prer$sample <- gsub("_\\d", "", prer$sample)

pres <- scRNA@meta.data %>% 
  rownames_to_column("cell") %>% 
  select(cell, drugSens, sample) %>% 
  filter(
    drugSens == "pre-sensitive" &
      !is.na(drugSens)
  )

pres$sample <- gsub("_\\d", "", pres$sample)

# Extract expression
exp <- subset(
  scRNA,
  cells = c(prer$cell, pres$cell))
exp <- exp@assays$RNA@data

samples <- c("NK", "Carbo", "Olaparib", "ctrl")
# Permutation
# cl <- getMPIcluster()
funtorun <- function(i, prer, pres, exp){
  samples <- c("NK", "Carbo", "Olaparib", "ctrl")
  prer_permu <- prer %>% 
    mutate(sample_permu = sample(prer$sample, replace = FALSE))
  
  pres_permu <- pres %>% 
    mutate(sample_permu = sample(pres$sample, replace = FALSE))
  
  classes_permu <- rbind.data.frame(prer_permu, pres_permu)
  
  fc <- NULL
  for (s in samples){
    pres_cell <- classes_permu$cell[
      which(classes_permu$drugSens == "pre-sensitive" & 
              classes_permu$sample_permu == s)
    ]
    pres_mean_exp <- log2(
      rowMeans(
        expm1(
          exp[, pres_cell]
        ),
        na.rm = TRUE
      ) + 1
    )
    prer_cell <- classes_permu$cell[
      which(classes_permu$drugSens == "pre-resistant" & 
              classes_permu$sample_permu == s)
    ]
    prer_mean_exp <- log2(
      rowMeans(
        expm1(
          exp[, prer_cell]
        ),
        na.rm = TRUE
      ) + 1
    )
    fc <- cbind(fc, prer_mean_exp - pres_mean_exp)
  }
  colnames(fc) <- paste0(samples, "_iter", i)
  fc <- as.data.frame(fc, stringsAsFactors = FALSE)
  fc$gene <- rownames(fc)
  return(fc)
}

fc_table <- lapply(
  1:p,
  function(i) funtorun(i=i, prer = prer, pres = pres, exp = exp)
)
# fc_table <- clusterApply(cl, 1:p, funtorun(prer = prer, pres = pres, exp = exp))
# stopCluster(cl)

fc_table <- reduce(
  fc_table,
  function(x, y){
    return(dplyr::full_join(x, y, by = "gene"))
  }
) 

fc_table <- fc_table %>% 
  left_join(
    scRNA@misc$gene_meta,
    by = c("gene" = "ensembl")
  ) %>% 
  select(gene, symbol, everything())
  

write.csv(fc_table, paste0(output_dir, "permutation_fc_table.csv"), row.names = FALSE)

# Observed value
fc_observe <- NULL
fc_observe_table <- NULL
classes_observe <- rbind.data.frame(prer, pres)
for (s in samples){
  pres_cell <- classes_observe$cell[
    which(classes_observe$drugSens == "pre-sensitive" &
            classes_observe$sample == s)
  ]
  pres_mean_exp <- log2(
    rowMeans(
      expm1(
        exp[, pres_cell]
      ),
      na.rm = TRUE
    ) + 1
  )
  prer_cell <- classes_observe$cell[
    which(classes_observe$drugSens == "pre-resistant" &
            classes_observe$sample == s)
  ]
  prer_mean_exp <- log2(
    rowMeans(
      expm1(
        exp[, prer_cell]
      ),
      na.rm = TRUE
    ) + 1
  )
  tmp_fc <- prer_mean_exp - pres_mean_exp
  if (is.null(fc_observe)){
    fc_observe <- tmp_fc
  } else {
    fc_observe <- cbind(fc_observe, tmp_fc)
  }
}

colnames(fc_observe) <- samples

fc_observe <- as.data.frame(fc_observe, stringsAsFactors = FALSE)
fc_observe <- fc_observe %>%
  rownames_to_column("gene") %>%
  left_join(
  scRNA@misc$gene_meta,
  by = c("gene" = "ensembl")
)

write.csv(fc_observe, paste0(output_dir, "permutation_fc_table_observation.csv"), row.names = FALSE)
fc_observe <- read.csv(paste0(output_dir, "permutation_fc_table_observation.csv"), stringsAsFactors = FALSE)

# p-value
combo <- combn(samples, 2)

p_value <- NULL
for (i in 1:ncol(combo)){
  # confidence interval of fc1-fc2
  tmp <- fc_table[, startsWith(colnames(fc_table), combo[1, i])] -
    fc_table[, startsWith(colnames(fc_table), combo[2, i])]
  tmp_observe <- fc_observe[, startsWith(colnames(fc_observe), combo[1, i])] -
    fc_observe[, startsWith(colnames(fc_observe), combo[2, i])]
  tmp_p <- rowSums(abs(tmp) > abs(tmp_observe))/p
  tmp_p <- data.frame(tmp_p, stringsAsFactors = FALSE)
  colnames(tmp_p) <- paste0(c(combo[, i], "p_value"), collapse = "_")
  tmp_p$gene <- fc_table$gene
  tmp_p$symbol <- fc_table$symbol
  if (is.null(p_value)){
    p_value<- tmp_p
  } else {
    p_value <- full_join(p_value, tmp_p, by = c("gene", "symbol"))
  }
}

p_value <- p_value %>% 
  select(gene, symbol, everything())

write.csv(p_value, paste0(output_dir, "permutation_p_value.csv"), row.names = FALSE)
