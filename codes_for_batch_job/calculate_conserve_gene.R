library(Seurat)
library(ggpubr)
library(tidyverse)
library(doParallel)
source("./auxiliary_codes.R")

merge <- readRDS("../data/output/KURAMOCHI_before_treatment/all_bt_sample_merge_filtered_gene_sisters.RDS")
output_dir <- "../data/output/KURAMOCHI_before_treatment/sister_conserve_gene/"
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# Extract twins from pooled pre-treatment samples

twins <- data.frame(cell_barcode = colnames(merge), 
                     sisters = merge$sisters_sample,
                     stringsAsFactors = FALSE)
twins <- twins[which(merge$true_sister), ]
twins <- twins[twins$sisters %in% names(table(twins$sisters)[table(twins$sisters) == 2]), ]
group <- unique(twins$sisters)
sisters <- NULL
for (g in group) {
  tmp <- data.frame(group = g,
                sister1 = twins$cell_barcode[twins$sisters == g][1],
                sister2 = twins$cell_barcode[twins$sisters == g][2],
                stringsAsFactors = FALSE)
  sisters <- rbind.data.frame(sisters, tmp)
}

cells_sister <- c(sisters$sister1, sisters$sister2)
n_pair <- nrow(sisters)
bt_twins <- subset(merge, cells = twins$cell_barcode)

# Extract gene expression table for twins

expression_mat <- as.matrix(GetAssayData(bt_twins, slot = "data"))
expression_mat <- log2(expm1(expression_mat)+1)
expression_mat <- expression_mat[, cells_sister]
expression_mat <- expression_mat[which(rowSums(expression_mat) != 0), ]
identical(colnames(expression_mat)[1:n_pair], sisters$sister1)
identical(colnames(expression_mat)[(n_pair + 1):ncol(expression_mat)], sisters$sister2)
dim(expression_mat)

log2FC_mat_sister <- abs(expression_mat[, 1:n_pair] - 
  expression_mat[, (n_pair + 1) : (2 * n_pair)])
cl <- makeCluster(40) # Specify the number of cores/clusters
registerDoParallel(cl)
t <- foreach(i = 1:nrow(expression_mat), .combine = "rbind.data.frame", .verbose = T) %dopar% {
  tmp <- NULL
  p <- NULL
  res <- NULL
  # Calculate the absolute difference of the gene i's expression on all the cell pairs 
  tmp <- as.matrix(dist(expression_mat[i,], method = "manhattan"))
  
  for (j in 1:nrow(sisters)){
    tmp[sisters$sister1[j], sisters$sister2[j]] <- NA
    tmp[sisters$sister2[j], sisters$sister1[j]] <- NA
  }
  tmp <- na.omit(tmp[upper.tri(tmp, diag = FALSE)])
  
  p <- t.test(log2FC_mat_sister[i, ], tmp)
  res <- data.frame(gene = rownames(expression_mat)[i],
                          mean_log2FC_sister = p$estimate["mean of x"],
                          mean_log2FC_non_sister = p$estimate["mean of y"],
                          log2FC_dif = p$estimate["mean of y"] - p$estimate["mean of x"],
                          t_statistic = p$statistic,
                          p_value = p$p.value,
                   stringsAsFactors = FALSE)
  return(res)
}
stopCluster(cl)

t$adj_p_value <- p.adjust(t$p_value, method = "fdr")
t <- dplyr::right_join(as.data.frame(merge@misc$gene_meta), t, by = c("ensembl" = "gene"))

write.csv(
  t, 
  paste0(output_dir, "/sister_conserve_gene_all_bt_sample_twin.csv"), 
  row.names = FALSE
)

# Conserved gene expression distribution

sum_exp <- read.csv(
  "../data/output/KURAMOCHI_before_treatment/genes_exp_sum_across_bt.csv",
  stringsAsFactors = FALSE)
sister_conserve_gene <- read.csv(
  paste0(output_dir, "sister_conserve_gene_all_bt_sample_twin.csv"),
  stringsAsFactors = FALSE
)

plot_table <- sister_conserve_gene %>% 
  mutate(conserve_gene = adj_p_value < 0.05) %>% 
  left_join(sum_exp, by = c("ensembl" = "gene"))
write.csv(plot_table, paste0(output_dir, "sum_exp_sister_conserve_gene_all_bt_sample_twin.csv"))
plot_table <-  plot_table %>%  
  select(ensembl, `Sister-conserved gene` = conserve_gene, sum_expression, mean_expression)

p <- ggplot(plot_table, aes(x = sum_expression, fill = `Sister-conserved gene`, color = `Sister-conserved gene`)) +
  geom_histogram(aes(y=..count..), position="identity", alpha=0.5, bins = 50) +
  scale_fill_manual(values = c("#CC331155", "#317DF755")) +
  scale_color_manual(values = c("#CC3311", "#317DF7")) +
  labs(x = "Sum of gene expression across before treatment samples",
       y = "Count of genes") +
  theme_classic()
p
ggsave(paste0(output_dir, "sister_conserve_gene_exp_distribution.pdf"), p)
ggsave(paste0(output_dir, "sister_conserve_gene_exp_distribution.png"),
       p, dpi = 300)
ggsave(paste0(output_dir, "sister_conserve_gene_exp_distribution.svg"), p)

# Extract genes in 2nd pin of the histogram.

plot_table_bin2 <- plot_table %>% 
  arrange(sum_expression)
group1 <- ggplot_build(p)$data[[1]] %>%
  filter(group == 1)
group2 <- ggplot_build(p)$data[[1]] %>%
  filter(group == 2)
n_pin_1 <- group1$y[1] + group2$y[1] # 618 + 19982 in old data
n_pin_2 <- group1$y[2] + group2$y[2] # 1742 + 1657 in old data
plot_table_bin2 <- plot_table_bin2[(n_pin_1 + 1):(n_pin_1 + n_pin_2), ]
plot_table_bin2 <- plot_table_bin2 %>% 
  left_join(select(sister_conserve_gene, ensembl, symbol), by = "ensembl")
write.csv(plot_table_bin2, "../data/output/KURAMOCHI_before_treatment/sister_conserve_gene/sum_exp_sister_conserve_gene_second_pin.csv",
          row.names = FALSE)
table(plot_table_bin2$`Sister-conserved gene`)

sessionInfo()
