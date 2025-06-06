#' Compute Geary's C for Multiple Samples
#' This function computes Geary's C score for multiple samples. It generates an upset plot to visualize the shared genes across samples.
#'
#' @param Expression_file A list containing expression data for each sample.
#' @param Meta_data A list containing metadata for each sample.
#' @param Output_path A string specifying the directory where output files will be saved.
#' @param Method A string indicating the method to be used for computation.
#' @param Tissue A string specifying the tissue type for the analysis.
#' @param pvalue A numeric value indicating the p-value threshold for selecting genes. Default is 0.01.
#'
#' @return A list containing the selected genes for each sample.
#' @import dplyr
#' @export
#'

Geary_C_score_multiple = function(Expression_file,Meta_data,Output_path,Method,Tissue,pvalue=0.01) {
  fractions = Fraction_multiple_samples(Expression_file,Meta_data)  
  Geary_list = list()
  samples = unique(Meta_data$sample)
  for (sample in samples) {
    Expression_file_subset = fractions$Expression_files[[sample]]
    Meta_data_subset = fractions$Meta_data_files[[sample]]
    print(sample)
    Geary_computation = Geary_C_score(Expression_file_subset,Meta_data_subset,Output_path = paste0(Output_path,sample,'_'),Method,paste0(Tissue,'_',sample),pvalue=pvalue)
    Geary_list[[sample]] = Geary_computation$Selected_genes
  }
  
  Shared_genes = unique(unlist(Geary_list))
  #Upset plot
  score_df = as.data.frame(Shared_genes)
  for (sample in samples) {
    score_df[[paste0('-',sample)]] <- as.integer(score_df$Shared_genes %in% Geary_list[[sample]])
  }
  pdf(paste(Output_path, Method, "_", Tissue, "_geary_samples_upset-plot.pdf", sep = ""), width = 6.5, height = 6.5, useDingbats = FALSE)
  p <- upset(score_df, sets = paste0('-',samples),
             order.by = "freq", empty.intersections = "off", text.scale = 1.5, matrix.color = "#56B4E9",
             main.bar.color = "grey", mainbar.y.label = 'Number of shared genes')
  print(p)
  dev.off()
  
  return(Geary_list)
}

