#' projection_visualization
#'
#' Visualize the results of projection step and inspect the major cell types in the samples and the number of shared cell types between the samples
#'
#'
#' @author Mengjie Chen
#' @param object A dmatch class object
#' @param filename The path and name of the output png file
#' @param TopCellLineNumber Keep only TopCellLineNumber primary cell lines which are highly correlated with any cell in the sampels. Regress the Pearson Correlation coefficients between the rest of primary cell lines with that specific cell in the samples to zero   
#' @param ShowCellNumber Keep only the primary cell lines which are highly correlated with more than ShowCellNumber cells in the samples 
#' @param dist.method The distance metric for calculating the distance between cells 
#' @param hclust.method The agglomeration method to be used for hierarchical cluster analysis
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information
#' @export
projection_visualization <- function(object, filename = NULL, TopCellLineNumber = 5, ShowCellNumber = 20, dist.method="euclidean", hclust.method="ward.D"){
  
  Projected <- object@Projection$cor.mat
  ReferenceNames <- object@Projection$ReferenceNames
  SampleNames <- object@Projection$SampleNames
  
  weights.mat <- apply(Projected, 1, function(x){
    y <- x
    cutoff <- x[order(x, decreasing=T)[TopCellLineNumber]]
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
  
  png(filename, res = 400, height = 8, width = 8, unit = "in")
  dist.method<-dist.method
  hclust.method<-hclust.method
  aa <- heatmap.2(kkk[flag, ], trace = "none", col = palette.gr.marray2, symbreaks = F,
                  labRow = ReferenceNames[flag], labCol = NA,  ColSideColors = colorlist[as.numeric(object@metadata$batch)],
                  key = F, margins = c(8, 15), distfun=function(x) dist(x,method = dist.method), hclustfun=function(x) hclust(x,method= hclust.method))
  dev.off()
  # bb <- rev(aa$colInd)
  
  bbb = kkk[flag, ]
  rownames(bbb) = ReferenceNames[flag]
  colnames(bbb) = SampleNames
  object@Projection.visualization<-list("Weight.max" = bbb, "ProjectionHeatmap" = aa)
  
  return(object)
  
}