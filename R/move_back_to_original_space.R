#' move_back_to_original_space
#'
#' Move cells in the PC space to the original high-dimensional space
#'
#'
#' @author Mengjie Chen
#' @param object A dmatch class object
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information
#' @export
move_back_to_original_space <- function(object){
  batch.id<-object@batch.id
  Final1 <- t(object@PCA$Raw[object@PCA$batch.id.forPC==batch.id[1],]) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Corrected - object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[1],])
  Final2 <- t(object@PCA$Raw[object@PCA$batch.id.forPC==batch.id[2],]) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Corrected - object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[2],])
  object@outcome.list<-list(Final1,Final2)
  return(object) 
}







