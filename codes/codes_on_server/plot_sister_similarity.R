library(ggplot2)
library(ggpubr)
source("./auxiliary_codes.R")
# without filter
p <- readRDS("../data/output/KURAMOCHI_before_treatment/sister_distance.RDS")
p1 <- PlotSisterSimilarity(distance = p$distance,
                          similarity = p$similarity,
                                 n_gene = p$n_gene,
                                 cutoff = NULL)
ggsave(paste0("../data/output/KURAMOCHI_before_treatment/", p$distance, "_", p$n_gene, "_genes.png"), plot = p1, device = "png", width = 6, height = 6, dpi = 300)
ggsave(paste0("../data/output/KURAMOCHI_before_treatment/", p$distance, "_", p$n_gene, "_genes.pdf"), plot = p1, device = "pdf", width = 6, height = 6)
# ggsave("../data/output/KURAMOCHI_before_treatment/euclidean_1000_genes_filter_sisters.png", plot = p1, device = "png", dpi = 600)
