### Curve

```{r}
sub_output_dir <- "./data/output/Figure4/curves/"
if (!dir.exists(sub_output_dir)) {
  dir.create(sub_output_dir)
}
# First batch
batches <- c(
  "synergyfinder_object_with_score_first_batch.RDS",
  "synergyfinder_object_with_score_second_batch.RDS",
  "synergyfinder_object_with_score_third_batch.RDS"
)
normalized <- TRUE

for (b in batches){
  data <- readRDS(paste0(combo_sub_dir, b))
  data$drug_pairs$drug1[grepl("NK", data$drug_pairs$drug1)] <- "NK"
  synergys <- c("HSA", "ZIP", "Loewe", "Bliss")
  
  for (plot_block in unique(data$drug_pairs$block_id)){
    drug_anno <- data$drug_pairs[data$drug_pairs$block_id == plot_block, ] %>% 
      dplyr::select(drug1, drug2, conc_unit1, conc_unit2)
    
    response <- data$response %>% 
      dplyr::filter(block_id == plot_block)
    
    drug2_concs <- unique(response$conc2)
    drug2_concs <- drug2_concs[order(drug2_concs)]
    # drug1_response <- list()
    models <- list()
    for (c in drug2_concs) {
      tmp <- response %>% 
        dplyr::filter(conc2 == c) %>% 
        dplyr::select(dose = conc1, response) %>% 
        dplyr::arrange(dose)
      if (normalized & c != "0") {
        control <- response %>% 
          dplyr::filter(conc2 == "0" & conc1 == "0") %>% 
          dplyr::select(response) %>% 
          unlist() %>% 
          mean()
        tmp$response <-  tmp$response - control
      }
      # drug1_response[[as.character(c)]] <- tmp
      tmp_model <- suppressWarnings(
        FitDoseResponse(
          tmp,
          Emin = NA,
          Emax = NA
        )
      )
      models[[as.character(c)]] <- tmp_model
    }
    
    plot_title <- paste(
      drug_anno$drug1, " and ", drug_anno$drug2
    )
    
    # plot the curve for the drug
    # For all of R's graphical devices, the default text size is 12 points but it
    # can be reset by including a pointsize argument to the function that opens
    #the graphical device. From ?pdf:
    # curve_pred <- data.frame(
    #   dose = seq(0, max(single_drug_data$dose), length.out = 700),
    #   response = PredictModelSpecify(
    #     models[[i]], 
    #     seq(0, max(single_drug_data$dose), length.out = 700)
    #   ),
    #   stringsAsFactors = FALSE
    # )
    # p1 <- ggplot(single_drug_data) +
    #   geom_point(aes(log10(dose), response)) +
    #   geom_line(aes(log10(dose), response), data = curve_pred) + 
    #   scale_y_continuous(trans = log10_trans())
    # p1
    colors <- c("#cacaca", "#979797", "#656565", "#323232", "#000000")
    names(colors) <- as.character(drug2_concs)
    ylim <- range(response$response)
    legend_position <- data.frame(
      x = 10 ^ -100,
      y = ylim[2] - seq(from = 5, length.out = length(models), by = 5)
    )
    rownames(legend_position) <- names(models)
    if (normalized) {
      pdf(file = paste0(sub_output_dir, drug_anno$drug1, "_", drug_anno$drug2, "_normalized_dose_response_curve.pdf"))
    } else {
      pdf(file = paste0(sub_output_dir, drug_anno$drug1, "_", drug_anno$drug2, "_dose_response_curve.pdf"))
    }
    
    for (i in names(models)){
      model_type <- FindModelType(models[[i]])
      if (model_type == "LL.4"){
        # Set break point of x-axis (numbers smaller than it are set as 0)
        bp <- round(min(log10(models[[i]]$origData$dose[models[[i]]$origData$dose > 1e-10])))-1
        max <- round(max(log10(models[[i]]$origData$dose)))
        step <- (max - bp)%/%4
        if (step < 1){
          step <- 1
        }
        xt <- 10 ^ seq(bp, max, by = step)
        xtlab <- xt
        xtlab[1] <- 0
      } else { # model type is L.4
        # Set break point of x-axis (numbers smaller than it are set as 0)
        bp <- round(min(log10(models[[i]]$origData$dose[models[[i]]$origData$dose > 1e-10])))-1
        # if (bp == 0) {
        #   bp <- -1
        # }
        max <- round(max(log10(models[[i]]$origData$dose)))
        step <- (max - bp)%/%4
        if (step < 1){
          step <- 1
        }
        # Set x-axis tick markders
        xt <- 10 ^ seq(bp, max, by = step)
        xtlab <- xt
        xtlab[1] <- 0
        models[[i]]$dataList$dose <- exp(models[[i]]$dataList$dose)
        models[[i]]$dataList$dose[models[[i]]$dataList$dose < 10^bp] <- 10^bp
      }
      
      # if (is.null(grid)){
      #   grid_exp <- NULL
      # } else if (grid){
      #   grid_exp <- expression(
      #     {
      #       graphics::grid(nx = NA, ny = NULL, col = "#DFDFDF", lty = 1)
      #       graphics::abline(col = "#DFDFDF", v = xt)
      #     }
      #   )
      # } else {
      #   grid_exp <- NULL
      # }
      grid_exp <- NULL
      # if (plot_new) {
      #   graphics::plot.new()
      #   grDevices::dev.control("enable")
      # }
      
      # suppressWarnings(graphics::par(plot_setting))
      
      # Plot dots
      # graphics::plot(
      #   x = models[[i]],
      #   xlab = paste0("Concentration (", drug_anno$unit, ")"),
      #   ylab = "Inhibition (%)",
      #   type = "obs",
      #   log = "x",
      #   col = "red",
      #   pch = 16,
      #   cex.axis = 1,
      #   panel.first = eval(grid_exp),
      #   xttrim = FALSE,
      #   conName = NULL,
      #   bp = 10 ^ bp,
      #   xt = xt,
      #   xtlab = xtlab
      # )
      
      # Plot curve
      if (i == "0"){
        graphics::plot(
          x = models[[i]],
          xlab = paste0(drug_anno$drug1, " (", drug_anno$conc_unit1, ")"),
          ylab = "Inhibition (%)",
          type = "none",
          log = "x",
          col = colors[[i]],
          cex.axis = 1,
          add = FALSE,
          lwd = 3,
          xttrim = FALSE,
          conName = NULL,
          bp = 10 ^ bp,
          xt = xt,
          xtlab = xtlab,
          ylim = ylim#,
          # legend = TRUE,
          # legendText = paste0(drug_anno$drug2, ": ", i, " ", drug_anno$conc_unit2),
          # legendPos = unlist(legend_position[i, ])
        )
      } else {
        graphics::plot(
          x = models[[i]],
          type = "none",
          log = "x",
          col = colors[[i]],
          cex.axis = 1,
          add = TRUE,
          lwd = 3,
          xttrim = FALSE,
          conName = NULL,
          bp = 10 ^ bp,
          xt = xt,
          xtlab = xtlab#,
          # legend = TRUE,
          # legendText = paste0(drug_anno$drug2, ": ", i, " ", drug_anno$conc_unit2),
          # legendPos = unlist(legend_position[i, ])
        )
      }
      
      # Plot title
      
      # graphics::mtext(
      #   plot_subtitle,
      #   cex = 7/9 * graphics::par()$cex.main * text_size_scale
      # )
      
      # if (record_plot) {
      #   p <- grDevices::recordPlot()
      #   grDevices::dev.off()
      #   print(p)
      #   return(p)
      # } else {
      #   return(NULL)
      # }
    }
    legend(
      "topleft", 
      legend=paste0(names(models), " ", drug_anno$conc_unit2),
      pch=15,
      fill = colors,
      title=drug_anno$drug2,
      pt.cex = 0.2
    )
    graphics::title(plot_title)
    grDevices::dev.off()
  }
}


```
