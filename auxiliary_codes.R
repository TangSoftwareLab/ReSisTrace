PlotQC <- function(object){
  
}
CalculateSisterSimilarity <- function(expression_mat, 
                                      n_gene = NULL, 
                                      sister_table, similarity = "euclidean"){
  if (!is.null(n_gene)){
    variance <- apply(expression_mat, 1, var)
    variance <- sort(variance, decreasing = TRUE)
    expression_mat <- expression_mat[which(rownames(expression_mat) %in% names(variance)[1:n_gene]), ]
  } else {
    n_gene = nrow(expression_mat)
  }
  
  sister_d <- rep(NA, nrow(sister_table))
  names(sister_d) <- sister_table$group
  if (similarity == "euclidean"){
    d <- wordspace::dist.matrix(t(expression_mat), method = "euclidean")
    d <- as.matrix(d)
  } else if (similarity == "pearson"){
    d <- cor(expression_mat, method = "pearson")
  } else if (similarity == "spearman"){
    d <- cor(expression_mat, method = "spearman")
  } else {
    stop("The similarity value ", similarity, " is not available! Available ",
         "values are: euclidean, pearson, and spearman.")
  }
  
  for (i in 1:nrow(sister_table)){
    sister_d[i] <- d[sister_table$sister1[i], sister_table$sister2[i]]
    d[sister_table$sister1[i], sister_table$sister2[i]] <- NA
    d[sister_table$sister2[i], sister_table$sister1[i]] <- NA
  }
  
  distance <- data.frame(similarity = c(na.omit(d[upper.tri(d, diag = FALSE)]), 
                                        sister_d),
                         stringsAsFactors = FALSE)
  distance$group <- c(rep("non-sisters", nrow(distance)-nrow(sister_table)), 
                      rep("sisters", nrow(sister_table)))
  res <- list(
    distance = distance,
    n_gene = n_gene,
    similarity = similarity,
    sister_distance = sister_d)
  return(res)
}

PlotSisterSimilarity <- function(distance,
                                 similarity = "euclidean",
                                 n_gene,
                                 cutoff = NULL){
  # Violin plot
  means <- aggregate(similarity ~  group, distance, mean)
  m <- switch (similarity,
               "euclidean" = "Euclidean Distance",
               "pearson" = "Pearson Coefficient",
               "spearman" = "Spearman Coefficient"
  )
  distance$group[distance$group == "non-sisters"] <- "Random cell pair"
  distance$group[distance$group == "sisters"] <- "Sister cell pair"
  p <- ggplot(data = distance, aes(x = group, y = similarity, fill = group)) +
    geom_violin(trim = FALSE) + 
    labs(title=paste0(m, " (genes: ", n_gene, ")"),x="", y = m) +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
  if (!is.null(cutoff)) {
    p <- p  + 
      geom_segment(aes(y = cutoff, x = 1.75, xend = 2.25, yend = cutoff))
  }
  p <- p + 
    scale_fill_manual(values=c("#317DF7", "#CC3311")) +
    theme_classic() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 20, face = "bold"),
          text = element_text(size = 20)) + 
    ggpubr::stat_compare_means(method = "t.test", vjust = 2, label.x.npc = "center") +
    geom_text(data = means, aes(label = signif(similarity, 4), y = similarity), 
              nudge_x = 0.15)
  return(p)
}


#' Calculate congruency level
#' 
#' @description This function calculates the 'congruncy level' from departure
#' cell to the destination cell. This score measures the similarity the lineage
#' labels detected in two cells.
#' 
#' @details Congruency is defined as the number of lineage labels of a departuring cell to the destination cell. For example:
#' \itemize{
#'     \item cell a has lineage label: L1 + L2 + L3 + L4
#'     \item cell b has lineage label: L1 + L2 + L3 + L4 + L5
#'     \item cell c has lineage label: L1
#'     \item cell d has lineage label: L1 + L2 + L3 + L5
#' }
#' 
#' Here we get:
#' \itemize{
#'     \item CL~ab~ (Congruent Level from cell a to cell b) = 4/5 = 80%
#'     \item CL~cb~ (Congruent level from cell c to cell b) = 1/5 = 20%
#'     \item CL~ba~ (Congruent level from cell b to cell a) = 5/4 = 125%
#' }  
#' 
#' We define the congruent level of cells contain conflict label as 0%. For example:
#'   
#' CL~da~ (Congruent level from cell d to cell a) = 0%
#' 
#' @param cell_from a vector of characters. It contains the lineage labels'
#' sequences or IDs detected from the departure cell.
#' @param cell_to a vector of characters. It contains the lineage labels'
#' sequences or IDs detected from the destination cell.
#' @return a numeric value. It the congruency level of the inputed two cells (from
#' departure cell to destination cell).
#' @export
#'
#' @examples
#' 
CongruencyLevel <- function(cell_from, cell_to){
  if (all(cell_from %in% cell_to) | all(cell_to %in% cell_from)){
    cl <- length(cell_from)/length(cell_to) * 100
  } else {
    cl <- 0
  }
  return(cl)
}

# expression_mat: the gene expression matrix in log2(normalized_count + 1) scale containing gene names as rows and cell barcodes as columns.
# n_gene: number of top variable genes used to calculate the similarity.
# lineage_table: a data frame containing the group number of families and corresponding preR or R cell's barcodes.

CalculatePreR_RSimilarity <- function(expression_mat,
                                      n_gene = NULL, 
                                      families,
                                      similarity = "euclidean"){
  if (!is.null(n_gene)){
    variance <- apply(expression_mat, 1, var)
    variance <- sort(variance, decreasing = TRUE)
    expression_mat <- expression_mat[which(rownames(expression_mat) %in% names(variance)[1:n_gene]), ]
  } else {
    n_gene = nrow(expression_mat)
  }
  
  bt_mat <- expression_mat[, startsWith(colnames(expression_mat), "BT_")]
  at_mat <- expression_mat[, startsWith(colnames(expression_mat), "AT_")]
  if (similarity == "euclidean"){
    d <- wordspace::dist.matrix(t(bt_mat), t(at_mat), method = "euclidean")
    d <- as.matrix(d)
  } else if (similarity == "pearson"){
    d <- cor(t(bt_mat), t(at_mat), method = "pearson")
  } else if (similarity == "spearman"){
    d <- cor(t(bt_mat), t(at_mat), method = "spearman")
  } else {
    stop("The similarity value ", similarity, " is not available! Available ",
         "values are: euclidean, pearson, and spearman.")
  }
  
  family_d <- rep(NA, nrow(families))
  names(family_d) <- families$group
  
  for (i in 1:nrow(families)){
    family_d[i] <- d[families$BT[i], families$AT[i]]
    d[families$BT[i], families$AT[i]] <- NA
  }
  
  distance <- data.frame(
    similarity = c(na.omit(c(d)), family_d),
    stringsAsFactors = FALSE
  )
  
  distance$group <- c(
    rep("Random BT-AT pair", nrow(distance)-nrow(families)),
    rep("PreR-R pair", nrow(families))
  )
  
  res <- list(
    distance = distance,
    n_gene = n_gene,
    similarity = similarity,
    family_distance = family_d)
  return(res)
}

PlotPreR_RSimilarity <- function(distance,
                                 similarity = "euclidean",
                                 n_gene,
                                 cutoff = NULL){
  # Violin plot
  means <- aggregate(similarity ~  group, distance, mean)
  m <- switch (similarity,
               "euclidean" = "Euclidean Distance",
               "pearson" = "Pearson Coefficient",
               "spearman" = "Spearman Coefficient"
  )
  p <- ggplot(data = distance, aes(x = group, y = similarity, fill = group)) +
    geom_violin(trim = FALSE) + 
    labs(title=paste0(m, " (genes: ", n_gene, ")"),x="", y = m) +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
  if (!is.null(cutoff)) {
    p <- p  + 
      geom_segment(aes(y = cutoff, x = 1.75, xend = 2.25, yend = cutoff))
  }
  p <- p + 
    scale_fill_manual(values=c("#317DF7", "#CC3311")) +
    theme_classic() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, face = "bold")) + 
    ggpubr::stat_compare_means(method = "t.test", vjust = 2, label.x.npc = "center") +
    ggpubr::stat_compare_means(method = "wilcox.test", vjust = 4, label.x.npc = "center") +
    geom_text(data = means, aes(label = signif(similarity, 4), y = similarity), 
              nudge_x = 0.15)
  return(p)
}
