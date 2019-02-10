#' cut_groups
#'
#' Cut cells in the samples into major cell types
#'
#'
#' @author Mengjie Chen
#' @param object A dmatch class object
#' @param K Number of cell types to cut
#' @param method The agglomeration method to be used for hierarchical cluster analysis   
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information
#' @export
cut_groups <- function(object, K, method="ward.D"){
  
  WeightMat <- object@Projection.visualization$Weight.max
  dd <- dist(t(WeightMat))
  ee <- hclust(dd, method)
  groups <- cutree(ee, k=K)
  groups<-as.data.frame(groups)
  rownames(groups)<-colnames(object@raw.data)
  colnames(groups)<-"CellType"
  object@metadata<-cbind(object@metadata,groups) 
  return(object)
}