#' move_back_to_original_space
#'
#' Move cells in the PC space to the original high-dimensional space, i.e., return the pairwise samples which are batch-effect corrected, centered, library-size normalized, and log transformed. 
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, outcome.list slot stores information for the corrected sample and the reference sample in the original high dimensional space
#' @export
move_back_to_original_space <- function(object){
  
  PCA<-object@PCA
  batch.id<-object@batch.id
  batch.id.forPC<-PCA$batch.id.forPC
  Data1 <- PCA$PCs[batch.id.forPC == batch.id[1],]
  Data2 <- PCA$PCs[batch.id.forPC == batch.id[2],]
  
  Raw1 <- object@PCA$Raw[batch.id.forPC == batch.id[1],]
  Raw2 <- object@PCA$Raw[batch.id.forPC == batch.id[2],]

  Final1 <- t(Raw1) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Corrected - Data1)
  Final2 <- t(Raw2) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Reference - Data2)
  object@outcome.list<-list("corrected_batch"=Final1,"reference_batch"=Final2)
  return(object) 
}







