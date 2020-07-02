#' cut_groups
#'
#' This step follows projection_visualization. After hclust() to cluster data and heatmap.2() to draw the heatmap in the projection_visualization step, this step uses cutree() to get subclusters, i.e., the major cell types in the samples. 
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object
#' @param K Number of cell types to cut
#' @param method The agglomeration method used by hclust() in the projection_visualization step for clustering data
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, cut_group slot stores information of celltypes for cells in the pairwise samples ("CellType")
#' @export
cut_groups <- function(object, K, method="ward.D"){
  
  WeightMat <- object@Projection.visualization$Weight.max
  dd <- dist(t(WeightMat))
  ee <- hclust(dd, method)
  groups <- cutree(ee, k=K)
  #cutoffs <- object@Projection.visualization$cutoffs
  object@cut_groups<-list("CellType"=groups, "batch.id.cut_groups"=object@Projection$batch.id.Proj.Vis)
  return(object)
}