#' Compute Variogram for Spatial Gene Expression Data
#'
#' This function computes the variogram for spatial gene expression data, allowing for
#' automatic padding, parallel processing, and optional saving of plots. It can handle
#' both data frames and summary expressions.
#'
#' @param Expression_file A data frame or matrix containing gene expression data.
#' @param Meta_data A data frame containing metadata for the cells, including their spatial coordinates.
#' @param automatic_pad A logical indicating whether to automatically determine the padding size. Default is TRUE.
#' @param n_pad An integer specifying the padding size. Default is 20. This is ignored if automatic_pad is TRUE.
#' @param save_plot A logical indicating whether to save the generated plots. Default is FALSE.
#' @param Output_path A character string specifying the path where output files (including plots) will be saved.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is 10.
#'
#' @return A list containing:
#'    - Selected_genes: A vector of selected genes based on spatial variability.
#'    - Parameters: A data frame with parameters of the variogram model for each gene, including:
#'        - Sill
#'        - Range
#'        - Nugget
#'        - Alpha
#'        - Total variance
#'        - R-squared
#'        - C/C+n (ratio of explained variance) 
#'        - AIC
#'        - Model
#'        - Variogram values
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom data.table as.data.table
#' @export
#' 



Compute_variogram = function(Expression_file,Meta_data,automatic_pad=TRUE,n_pad=20,save_plot=FALSE,Output_path,ncores= 10) {
  registerDoParallel(cores=ncores)
  if (automatic_pad) {
    n_pad = pad_definitor(Meta_data)
    cat('Selected automatic pad:', n_pad, '\n')
  }
  breaks_x <- seq(min(Meta_data$cell_centroid_x), max(Meta_data$cell_centroid_x) + n_pad, by = n_pad)
  breaks_y <- seq(min(Meta_data$cell_centroid_y), max(Meta_data$cell_centroid_y) + n_pad, by = n_pad)
  if (is.data.frame(Expression_file)) {
  
  # Process the data
  pos_exp <- Meta_data %>%
    dplyr::select(cell_centroid_x, cell_centroid_y) %>%
    bind_cols(Expression_file) %>%
    arrange(cell_centroid_x) %>%
    mutate(
      portion_x = cut(cell_centroid_x, breaks = breaks_x, right = FALSE),
      portion_y = cut(cell_centroid_y, breaks = breaks_y, right = FALSE)
    ) %>%
    mutate(across(-c(cell_centroid_x, cell_centroid_y, portion_x, portion_y), as.numeric))
  pos_exp_dt <- as.data.table(pos_exp)
  sum_by_cut <- pos_exp_dt[, lapply(.SD, sum, na.rm = TRUE), by = .(portion_x, portion_y), .SDcols = !c("portion_x", "portion_y")]
  colnames(sum_by_cut) <- gsub("/", "-", colnames(sum_by_cut))

  gene_names <- colnames(sum_by_cut[,-c(1,2,3,4)])
  }
  else {
    Expression_tibble <- summary(Expression_file) %>%
      as_tibble(.name_repair = "minimal")
    
    colnames(Expression_tibble) <- c("cell", "gene", "expression")
    Expression_tibble <- Expression_tibble %>%
      mutate(cell = cell, 
             gene = colnames(Expression_file)[gene]) 
    
    Meta_data_tibble <- Meta_data %>%
      mutate(cell = row_number()) %>% 
      dplyr::select(cell, cell_centroid_x, cell_centroid_y)
    gc()
    Expression_tibble <- Expression_tibble %>%
      left_join(Meta_data_tibble, by = "cell")
    
    Expression_tibble <- Expression_tibble %>%
      mutate(
        portion_x = cut(cell_centroid_x, breaks = breaks_x, right = FALSE),
        portion_y = cut(cell_centroid_y, breaks = breaks_y, right = FALSE)
      )
    gc()
    sum_by_cut <- Expression_tibble %>%
      group_by(portion_x, portion_y, gene) %>%
      summarise(expression = sum(expression, na.rm = TRUE), .groups = "drop")
    gc()
    sum_by_cut <- sum_by_cut %>%
      pivot_wider(names_from = gene, values_from = expression, values_fill = 0)
    colnames(sum_by_cut) <- gsub("/", "-", colnames(sum_by_cut))
    gene_names <- colnames(sum_by_cut[,-c(1,2)])
    
  }
  gc()
  ptime <- system.time({
    results <- foreach(gene = gene_names, .combine = 'rbind') %dopar% {
      process_gene(gene, sum_by_cut, n_pad, save_plot = save_plot, Output_path = Output_path)
    }
  })[3]
  ptime
  stopImplicitCluster()
  gc()
  #Extract spatially variable genes
  var_parameters <- as.data.frame(results)
  colnames(var_parameters) <- c('Sill', 'Range', 'Nugget(n)','Alpha','Total variance (C+n)', 'R2', 'C/C+n','AIC','Model','Vario_values') #
  var_parameters$gene_names <- gene_names
  rownames(var_parameters) <- gene_names
  var_parameters[, 1:8] <- lapply(var_parameters[, 1:8], as.numeric)
  var_parameters$Range_corrected <- ifelse(var_parameters$Model == 'Exponential',
                                           var_parameters$Range * log(2), 
                                           ifelse(var_parameters$Model == 'Finetuned exponential',
                                                  (var_parameters$Range * log(2))^(1/var_parameters$Alpha),
                                                     var_parameters$Range))
  
  gc()
  pdf(paste0(Output_path,"models_distribution.pdf"), useDingbats = FALSE, width = 8, height = 8)
  par(las=1, bty="l", mar=c(5, 4, 4, 2) + 0.1)
  
  all_models <- c('Constant','Exponential','Finetuned exponential')
  model_counts <- table(unlist(var_parameters$Model))
  model_percentages <- model_counts / nrow(var_parameters) * 100
  model_percentages <- as.data.frame(model_percentages)
  colnames(model_percentages) <- c("Model", "Percentage")
  complete_model_percentages <- data.frame(Model = all_models, Percentage = 0)
  complete_model_percentages <- bind_rows(model_percentages,anti_join(complete_model_percentages,model_percentages,by='Model'))
  barplot(complete_model_percentages$Percentage, beside=TRUE,ylim = c(0, 100), col = 'red', main = "Fitted Models Distribution",
          xlab = "Models", ylab = "Percentage of fitted genes (%)",names.arg = complete_model_percentages$Model) 
  text(x=barplot(complete_model_percentages$Percentage, beside=TRUE, plot=FALSE),
       y=complete_model_percentages$Percentage,
       labels=round(complete_model_percentages$Percentage, 2),
       pos=3)
  Selected = gene_names[which(var_parameters$R2 > 0.85 & var_parameters$`C/C+n` > 0.3 & var_parameters$`Total variance (C+n)` > 0.5) ]  #
  dev.off()
  return(list(Selected_genes = Selected,Parameters = var_parameters))
}

