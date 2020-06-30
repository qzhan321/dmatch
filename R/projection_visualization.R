#' projection_visualization
#'
#' Visualize the results of projection step and inspect the major cell types in the samples and the number of shared cell types between the samples. Sparsity is introduced in the correlation matrix by only retaining the top 5 (default) primary cell lines (all others are set to zero). A minimum threshold for the correlations is provided via cor.threshold parameter. If the overall correlation between cells in the samples and the last TopCellLine in the reference panel is less than the threshold, a warning will be generated. Our method does not work well with those datasets.
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object.
#' @param filename The path and name of the output png file. Default is NULL, and a heatmap is printed on the screen directly.
#' @param TopCellLineNumber The number of primary cell lines to be kept which are highly correlated with any cell in the samples. Set the Pearson Correlation coefficients between the rest of primary cell lines with cells in the samples to zero.
#' @param cor.threshold If the correlation between some cells in the samples and the last TopCellLine in the reference is lower than the threshold, a warning will be generated.
#' @param ShowCellNumber Include only the primary cell lines which are highly correlated with more than this amount of cells in the samples.
#' @param dist.method The distance metric for calculating the distance between cells.
#' @param hclust.method The agglomeration method used by hclust() in the projection_visualization step for clustering data
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, projection.visualization slot stores information for the weight matrix between cells in the samples and primary cell lines in the reference with enforeced sparsity, and a heat map.
#' @export
projection_visualization <- function(object, filename = NULL, TopCellLineNumber = 5, cor.threshold = 0.6, ShowCellNumber = 20, dist.method="euclidean", hclust.method="ward.D"){
  
  Projected <- object@Projection$cor.mat
  ReferenceNames <- object@Projection$ReferenceNames
  SampleNames <- object@Projection$SampleNames
  
  cutoffs <- c()
  weights.mat <- apply(Projected, 1, function(x){
    y <- x
    cutoff <- x[order(x, decreasing=T)[TopCellLineNumber]]
    if (cutoff < cor.threshold) {
      stop(paste("The correlation between some cells in the data and the last TopCellLine is low which will result in an unrealible identification of clusters in the data. Recomend to use less TopCellLines."))
    }
    cutoffs <- c(cutoffs, cutoff)
    y[which(x<cutoff)] <- 0
    y
  })
  
  flag <- apply(weights.mat, 1, function(x){
    length(x[x!=0]) >= ShowCellNumber
  })
  
  kkk <- as.matrix(weights.mat)
  
  require(gplots)
  palette.gr.marray2 <- colorRampPalette(c("white", "red"))(16)
  colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki",
                 "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",
                 "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
  if (is.null(filename)) {
    dist.method<-dist.method
    hclust.method<-hclust.method
    # heatmap.2(kkk[flag, ], trace = "none", col = palette.gr.marray2, symbreaks = F,
    #           labRow = ReferenceNames[flag], labCol = NA,  ColSideColors = colorlist[as.numeric(object@metadata$batch)],
    #           key = F, margins = c(8, 15), distfun=function(x) dist(x,method = dist.method), hclustfun=function(x) hclust(x,method= hclust.method))
    aa <- heatmap.2(kkk[flag, ], trace = "none", col = palette.gr.marray2, symbreaks = F,
                    labRow = ReferenceNames[flag], labCol = NA,  ColSideColors = colorlist[as.numeric(object@Projection$batch.id.update)],
                    key = F, margins = c(8, 15), distfun=function(x) dist(x,method = dist.method), hclustfun=function(x) hclust(x,method= hclust.method))
  } else {
    png(filename, res = 400, height = 8, width = 8, unit = "in")
    dist.method<-dist.method
    hclust.method<-hclust.method
    aa <- heatmap.2(kkk[flag, ], trace = "none", col = palette.gr.marray2, symbreaks = F,
                    labRow = ReferenceNames[flag], labCol = NA,  ColSideColors = colorlist[as.numeric(object@Projection$batch.id.update)],
                    key = F, margins = c(8, 15), distfun=function(x) dist(x,method = dist.method), hclustfun=function(x) hclust(x,method= hclust.method))
    dev.off()
    # bb <- rev(aa$colInd)
  }
  
  bbb = kkk[flag, ]
  rownames(bbb) = ReferenceNames[flag]
  colnames(bbb) = SampleNames
  object@Projection.visualization<-list("Weight.max" = bbb, "ProjectionHeatmap" = aa, "cutoffs" = cutoffs)
  
  return(object)
}
