#' Save Comparison Plots for Clustering and Annotation
#'
#' This function generates and saves comparison plots for clustering results and cell annotations.
#' It creates KNN and Scimilarity plots in both tissue atlas and UMAP formats.
#'
#' @param Clustering A vector indicating the clustering assignments for each cell.
#' @param Annotation_cells_renamed A vector indicating the annotations for each cell.
#' @param Umap_result A UMAP object containing UMAP coordinates for each cell.
#' @param Output_path A character string specifying the directory path where the output PDF plot will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @param scaling_factor A numeric value for scaling the point sizes in the plots. Default is 2.
#'
#' @return  NULL This function does not return a value; it saves PNG files.
#' @export
#' 
Save_comparison_plots = function(Clustering,Annotation_cells_renamed,Umap_result,Output_path,Method,Tissue,scaling_factor=2) {
  ppp_temp = spatstat::ppp(x = Meta_data$cell_centroid_x,y = Meta_data$cell_centroid_y,
                 window = owin(xrange = range(Meta_data$cell_centroid_x),yrange = range(Meta_data$cell_centroid_y)))
  ratio <- table(Clustering, Annotation_cells_renamed)
  filtered_table <- as.data.frame(ratio)
  filtered_table <- filtered_table %>%
    group_by(Clustering) %>%
    dplyr::filter(Freq == max(Freq))
  filtered_table <- filtered_table %>%
    arrange(Annotation_cells_renamed, desc(Freq))
  
  #### KNN Clustering
  N_cluster = length(filtered_table$Clustering)
  optimal_palette = rainbow(N_cluster)
  png(paste0(Output_path,Method,"_",Tissue,"_KNN_atlas.png"),width = 4000,height = 4000)
  plot(ppp_temp[1,],main="",bty="n",cex = 0)
  for (k in 1:length(filtered_table$Clustering)) {
    temp_cell_type = filtered_table$Clustering[k]
    points(ppp_temp[Clustering==temp_cell_type,],cex=scaling_factor,col=optimal_palette[k])
  }
  dev.off()
  
  ### SCIMILARITY
  del = NULL
  for (k in 2:nrow(filtered_table)) {
    if (filtered_table$Annotation_cells_renamed[k] == filtered_table$Annotation_cells_renamed[k-1]) {
      del  = c(del, k)
    }
  }
  optimal_palette2 <- optimal_palette[-del]
  png(paste0(Output_path,Method,"_",Tissue,"_Scimilarity_atlas.png"),width = 4000,height = 4000)
  plot(ppp_temp[1,],main="",bty="n",cex = 0)
  for (k in 1:length(unique(filtered_table$Annotation_cells_renamed))) {
    temp_cell_type = unique(filtered_table$Annotation_cells_renamed)[k]
    points(ppp_temp[Annotation_cells_renamed==temp_cell_type,], cex=scaling_factor, col=optimal_palette2[k])
  }
  dev.off()
  
  ###KNN
  names(optimal_palette) <- filtered_table$Clustering
  png(paste0(Output_path,Method,"_",Tissue,"_KNN_umap.png"),width = 1000,height = 1000)
  plot(Umap_result,pch=21,bg=optimal_palette[Clustering],cex=1,lwd=0.2)
  legend("topleft", legend = filtered_table$Clustering, fill =optimal_palette, title = "KNN", lty = 1,cex = 1, box.lty = 0.4 , box.lwd = 0.4)
  dev.off()
  ###Scimilarity
  names(optimal_palette2) <- unique(filtered_table$Annotation_cells_renamed)
  png(paste0(Output_path,Method,"_",Tissue,"_Scimilarity_umap.png"),width = 1000,height = 1000)
  plot(Umap_result,pch=21,bg=optimal_palette2[Annotation_cells_renamed],cex=1,lwd=0.2)
  legend("topleft", legend = unique(filtered_table$Annotation_cells_renamed), fill = optimal_palette2, title = "Scimilarity", lty = 1, cex = 0.8, box.lty = 0.5 , box.lwd = 1)
  dev.off()
  return(invisible())
}

