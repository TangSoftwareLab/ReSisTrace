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

TwoDoseResponseCurve <- function(
    data,
    plot_block = 1,
    drug_index = 1,
    other_drug_dose = 2,
    adjusted = TRUE,
    Emin = NA,
    Emax = NA,
    grid = TRUE,
    point_color = "#C24B40",
    curve_color = "black",
    point_color_2 = "blue",
    curve_color_2 = "black",
    text_size_scale = 1,
    plot_title = NULL,
    plot_subtitle = NULL,
    plot_setting = list(
      cex.lab = 1 * text_size_scale,
      mgp = c(2, 0.5, 0),
      font.main = 2,
      font.lab = 1,
      cex.main = 14 / 12 * text_size_scale,
      bty = "l",
      lwd = 1.5
    ),
    plot_new = TRUE,
    record_plot = TRUE,
    ...) {
  
  # 1. Check the input data
  # Data structure of 'data'
  if (!is.list(data)) {
    stop("Input data is not in list format!")
  }
  if (!all(c("drug_pairs", "response") %in% names(data))) {
    stop("Input data should contain at least tow elements: 'drug_pairs' and 
         'response'. Please prepare your data with 'ReshapeData' function.")
  }
  # Parameter 'plot_block'
  if (length(drug_index) != 1) {
    stop("The length of 'plot_block' parameter is not 1. Please chosed only one
         block for visualization.")
  }
  # Parameter 'drug'
  concs <- grep("conc\\d", colnames(data$response), value = TRUE)
  if (length(drug_index) != 1) {
    stop("The length of 'drug' parameter is not 1. Please chosed only one
         drug in the block for visualization.")
  } else if (!drug_index %in% sub("conc", "", concs)) {
    stop("The input drug_index '", drug_index, "' is not found in input data. ",
         "Available indexes are '", 
         paste(sub("conc", "", concs), collapse = "', '"),
         "'.")
  }
  
  # Annotation data
  drug_anno <- data$drug_pairs[data$drug_pairs$block_id == plot_block, ] %>% 
    dplyr::select(drug = paste0("drug", drug_index),
                  unit = paste0("conc_unit", drug_index))
  
  # Extract single drug dose response
  if (!adjusted) {
    response <- data$response %>% 
      dplyr::select(-response) %>% 
      dplyr::rename(response = response_origin) %>% 
      dplyr::filter(block_id == plot_block)
  } else {
    response <- data$response %>% 
      dplyr::filter(block_id == plot_block)
  }
  single_drug_data <- ExtractSingleDrug(response)
  single_drug_data <- single_drug_data[[paste0("conc", drug_index)]]
  # Fit model for the row drug
  drug_model <- suppressWarnings(
    FitDoseResponse(
      single_drug_data,
      Emin = Emin,
      Emax = Emax
    )
  )
  
  if (is.null(plot_subtitle)) {
    plot_subtitle <- paste(
      drug_anno$drug,
      "in Block",
      plot_block
    )
  }
  if (is.null(plot_title)) {
    plot_title <- "Dose-Response Curve"
  }
  
  # plot the curve for the drug
  # For all of R's graphical devices, the default text size is 12 points but it
  # can be reset by including a pointsize argument to the function that opens
  #the graphical device. From ?pdf:
  # curve_pred <- data.frame(
  #   dose = seq(0, max(single_drug_data$dose), length.out = 700),
  #   response = PredictModelSpecify(
  #     drug_model, 
  #     seq(0, max(single_drug_data$dose), length.out = 700)
  #   ),
  #   stringsAsFactors = FALSE
  # )
  # p1 <- ggplot(single_drug_data) +
  #   geom_point(aes(log10(dose), response)) +
  #   geom_line(aes(log10(dose), response), data = curve_pred) + 
  #   scale_y_continuous(trans = log10_trans())
  # p1
  
  model_type <- FindModelType(drug_model)
  if (model_type == "LL.4"){
    # Set break point of x-axis (numbers smaller than it are set as 0)
    bp <- round(min(log10(drug_model$origData$dose[drug_model$origData$dose > 1e-10])))-1
    max <- round(max(log10(drug_model$origData$dose)))
    step <- (max - bp)%/%4
    if (step < 1){
      step <- 1
    }
    xt <- 10 ^ seq(bp, max, by = step)
    xtlab <- xt
    xtlab[1] <- 0
  } else { # model type is L.4
    # Set break point of x-axis (numbers smaller than it are set as 0)
    bp <- round(min(log10(drug_model$origData$dose[drug_model$origData$dose > 1e-10])))-1
    # if (bp == 0) {
    #   bp <- -1
    # }
    max <- round(max(log10(drug_model$origData$dose)))
    step <- (max - bp)%/%4
    if (step < 1){
      step <- 1
    }
    # Set x-axis tick markders
    xt <- 10 ^ seq(bp, max, by = step)
    xtlab <- xt
    xtlab[1] <- 0
    drug_model$dataList$dose <- exp(drug_model$dataList$dose)
    drug_model$dataList$dose[drug_model$dataList$dose < 10^bp] <- 10^bp
  }
  
  if (is.null(grid)){
    grid_exp <- NULL
  } else if (grid){
    grid_exp <- expression(
      {
        graphics::grid(nx = NA, ny = NULL, col = "#DFDFDF", lty = 1)
        graphics::abline(col = "#DFDFDF", v = xt)
      }
    )
  } else {
    grid_exp <- NULL
  }
  
  if (plot_new) {
    graphics::plot.new()
    grDevices::dev.control("enable")
  }
  
  suppressWarnings(graphics::par(plot_setting))
  
  # Plot dots
  graphics::plot(
    x = drug_model,
    xlab = paste0("Concentration (", drug_anno$unit, ")"),
    ylab = "Inhibition (%)",
    type = "obs",
    log = "x",
    col = point_color,
    pch = 16,
    cex.axis = 1 * text_size_scale,
    panel.first = eval(grid_exp),
    xttrim = FALSE,
    conName = NULL,
    bp = 10 ^ bp,
    xt = xt,
    xtlab = xtlab,
    ...
  )
  
  # Plot curve
  graphics::plot(
    x = drug_model,
    type = "none",
    log = "x",
    col = curve_color,
    cex.axis = 1 * text_size_scale,
    add = TRUE,
    lwd = 3,
    xttrim = FALSE,
    conName = NULL,
    bp = 10 ^ bp,
    xt = xt,
    xtlab = xtlab,
    ...
  )
  
  # Plot title
  graphics::title(plot_title)
  graphics::mtext(
    plot_subtitle,
    cex = 7/9 * graphics::par()$cex.main * text_size_scale
  )
  
  if (record_plot) {
    p <- grDevices::recordPlot()
    grDevices::dev.off()
    print(p)
    return(p)
  } else {
    return(NULL)
  }
}

#' Extract Single Drug Dose Response
#'
#' \code{ExtractSingleDrug} extracts the dose-response values of single drug
#'  from a drug combination dose-response matrix.
#'
#' @param response A data frame. It must contain the columns: "conc1", "conc2",
#' ..., for the concentration of the combined drugs and "response" for the
#' observed \%inhibition at certain combination.
#'
#' @return A list contains several data frames each of which contains 2 columns:
#'   \itemize{
#'     \item \strong{dose} The concertration of drug.
#'     \item \strong{response} The cell's response (inhibation rate) to
#'       corresponding drug concertration.
#' }
#'
#' @author
#' \itemize{
#'   \item Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'   \item Jing Tang \email{jing.tang@helsinki.fi}
#' }
#' 
#' @export
#'
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' response <- data$response[data$response$block_id == 1,
#'                           c("conc1", "conc2", "response")]
#' single <- ExtractSingleDrug(response)
ExtractSingleDrug <- function(response) {
  concs <- grep("conc\\d", colnames(response), value = TRUE)
  single_drug <- vector("list", length(concs))
  names(single_drug) <- concs
  conc_sum <- rowSums(response[, concs])
  for (conc in concs) {
    index <- which(response[, conc] == conc_sum)
    single_drug[[conc]] <- data.frame(
      dose = unlist(response[index, conc]),
      response = response[index, "response"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  return(single_drug)
}