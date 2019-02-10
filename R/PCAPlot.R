#' PCAPlot
#'
#' Visualize the result of SVD
#'
#'
#' @author Mengjie Chen
#' @param PCA The list returned by fastSVD
#' @param PCs.to.plot Which two PCs to visualize
#' @param filename The path and name of the PCA plot returned
#' @return A png file for the visualization of SVD results
#' @export
PCAPlot<-function(PCA, PCs.to.plot, filename) {
  i<-PCs.to.plot[1]
  j<-PCs.to.plot[2]
  PC1 <- PCA$PCs[, i]
  PC2 <- PCA$PCs[, j]
  batch.id.forPC<-PCA$batch.id.forPC
  png(filename, res = 400, height = 8, width = 8, unit = "in")
  plot(PC1, PC2, type = "n")
  text(PC1, PC2, batch.id.forPC, col = batch.id.forPC, cex = 0.8)
  dev.off()
}




