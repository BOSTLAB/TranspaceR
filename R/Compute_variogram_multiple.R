#' Compute Variogram for Multiple Samples
#'
#' This function computes the variogram for multiple samples based on the provided expression data and metadata.
#' It generates an upset plot to visualize the shared genes across samples.
#'
#' @param Expression_file A data frame or file path containing the expression data for the samples.
#' @param Meta_data A data frame containing metadata for the samples, including sample identifiers.
#' @param automatic_pad A logical value indicating whether to automatically pad the data. Default is TRUE.
#' @param n_pad An integer specifying the number of padding points to use if automatic padding is enabled. Default is 20.
#' @param save_plot A logical value indicating whether to save the generated plot. Default is FALSE.
#' @param Output_path A character string specifying the directory where output files (including plots) will be saved.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is 10.
#'
#' @return A list containing the selected genes for each sample based on the computed variograms.
#'                
#' @export
Compute_variogram_multiple = function(Expression_file,Meta_data,automatic_pad=TRUE,n_pad=20,save_plot=FALSE,Output_path,ncores= 10) {
  fractions = Fraction_multiple_samples(Expression_file,Meta_data)  
  Variogram_list = list()
  samples = unique(Meta_data$sample)
  for (sample in samples) {
    Expression_file_subset = fractions$Expression_files[[sample]]
    Meta_data_subset = fractions$Meta_data_files[[sample]]
    print(sample)
    Variogram_computation = Compute_variogram(Expression_file_subset,Meta_data_subset,
                            automatic_pad=FALSE,n_pad=100,save_plot=TRUE,ncores=5,paste0(Output_path,sample,'_'))
    Variogram_list[[sample]] = Variogram_computation$Selected_genes
  }
  
  Shared_genes = unique(unlist(Variogram_list))
  #Upset plot
  score_df = as.data.frame(Shared_genes)
  for (sample in samples) {
    score_df[[paste0('-',sample)]] <- as.integer(score_df$Shared_genes %in% Variogram_list[[sample]])
  }
  pdf(paste(Output_path, Method, "_", Tissue, "_variogram_samples_upset-plot.pdf", sep = ""), width = 6.5, height = 6.5, useDingbats = FALSE)
  p <- upset(score_df, sets = paste0('-',samples),
             order.by = "freq", empty.intersections = "off", text.scale = 1.5, matrix.color = "#56B4E9",
             main.bar.color = "grey", mainbar.y.label = 'Number of shared genes')
  print(p)
  dev.off()
  
  return(Variogram_list)
}
