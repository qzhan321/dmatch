#' PCAPlot
#'
#' Visualize the result of SVD
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param PCA The output from fastSVD, i.e., the list returned by fastSVD function
#' @param PCs.to.plot A vector indicates which two PCs to plot
#' @param batchs.to.plot A vector indicates which batches of PCs to plot
#' @param filename The path and name of the PCA plot returned, default value is NULL and the PCA plot is printed directly 
#' @return Graphs the output of SVD results, colored by their batch id
#' @export
PCAPlot<-function(PCA, PCs.to.plot, batchs.to.plot, filename=NULL) {
  i<-PCs.to.plot[1]
  j<-PCs.to.plot[2]
  
  batch.id.forPC<-PCA$batch.id.forPC
  batch.id<-batch.id.forPC[batch.id.forPC %in% batchs.to.plot]
  
  xx <- PCA$PCs[batch.id.forPC %in% batchs.to.plot, i]
  yy <- PCA$PCs[batch.id.forPC %in% batchs.to.plot, j]
  
  if (is.null(filename)) {
    plot(xx, yy, type = "n", xlab = paste0("PC", i), ylab = paste0("PC", j))
    text(xx, yy, batch.id, col = batch.id, cex = 0.8)
  } else {
    png(filename, res = 400, height = 8, width = 8, unit = "in")
    plot(xx, yy, type = "n", xlab = paste0("PC", i), ylab = paste0("PC", j))
    text(xx, yy, batch.id, col = batch.id, cex = 0.8)
    dev.off()
  }
}




