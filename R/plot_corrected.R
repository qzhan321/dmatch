#' plot_corrected
#'
#' PC plots for corrected samples. 
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object.
#' @param filename The path and name of the png file. Default is NULL and the PC plot is printed directly.
#' @param pcs.plot Which two PCs to plot.
#' @param celltypes.plot Which cell clusters to plot.
#' @return Return graphs of samples in PC space after batch effects correction.
#' @export
plot_corrected <- function(object, filename = NULL, pcs.plot, celltypes.plot){
  
  SharedType <- object@run_alignment_by_2D.results$RefTypes
  
  Data2 <- object@run_alignment_by_2D.results$Reference
  Data1 <- object@run_alignment_by_2D.results$Original
  
  Labels1 <- object@run_alignment_by_2D.results$Labels1
  Labels2 <- object@run_alignment_by_2D.results$Labels2
  
  Corrected <- object@run_alignment_by_2D.results$Corrected
  
  pc.num1<-pcs.plot[1]
  pc.num2<-pcs.plot[2]
  
  if (is.null(filename)) {
    par(mfrow = c(1, 2))
    xlim <- c(min(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))-1, max(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))+1)
    ylim <- c(min(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))-1, max(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))+1)
    
    plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "Before Correction")
    for(i in 1:length(celltypes.plot)){
      X <- Data1[Labels1 == celltypes.plot[i], pcs.plot]
      x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
      points(x1, col = "black", pch = i)
      points(X, col = "red", pch = i)
    }
    
    plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "After Correction")
    for(i in 1:length(celltypes.plot)){
      x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
      X <- Corrected[Labels1 == celltypes.plot[i], pcs.plot]
      points(x1, col = "black", pch = i)
      points(X, col = "green", pch = i)
    }
    
  } else {
    png(filename, res = 400, height = 6, width = 10, unit = "in")
    par(mfrow = c(1, 2))
    xlim <- c(min(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))-1, max(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))+1)
    ylim <- c(min(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))-1, max(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))+1)
    
    plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "Before Correction")
    for(i in 1:length(celltypes.plot)){
      X <- Data1[Labels1 == celltypes.plot[i], pcs.plot]
      x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
      points(x1, col = "black", pch = i)
      points(X, col = "red", pch = i)
    }
    
    plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "After Correction")
    for(i in 1:length(celltypes.plot)){
      x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
      X <- Corrected[Labels1 == celltypes.plot[i], pcs.plot]
      points(x1, col = "black", pch = i)
      points(X, col = "green", pch = i)
    }
    dev.off()
  }
}


