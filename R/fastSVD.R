#' fastSVD
#'
#' Perform SVD
#'
#'
#' @author Mengjie Chen
#' @param sample.list A list of samples
#' @param nPC Number of PCs to keep
#' @param min.cells Filtering out the genes which are expressed in less than min.cells cells
#' @param min.genes Filtering out the cells which have less than min.genes genes expressed
#' @return A list consists of: PCs, Loadings, Centers, batch.id for samples while performing SVD, Raw representing the input (before centering) for SVD, cells used for SVD, genes used for SVD
#' @export

fastSVD <- function(samples.list, nPC = 30, min.cells=NULL, min.genes=NULL){
  alldata <- samples.list[[1]]
  batch.id <- rep(1,ncol(alldata))
  for (i in 2:length(samples.list)) {
    sample <- samples.list[[i]]
    alldata <- cbind(alldata, sample)
    batch.id <- c(batch.id, rep(i, ncol(sample)))
  }  
  if (is.null(min.cells)) {
    alldata<-alldata
  }  else {
    alldata<-alldata[rowSums(alldata>0)>=min.cells,]
  }
  if (is.null(min.genes)) {
    alldata<-alldata
  } else {
    batch.id<-batch.id[colSums(alldata)>=min.genes]
    alldata<-alldata[,colSums(alldata)>=min.genes]
  }
  
  alldata_norm <- apply(alldata, 2, function(x) x/sum(x))
  alldata_libsize <- apply(alldata_norm, 2, function(x) x*1000000)
  alldata_log <- apply(alldata_libsize, 2, function(x) log(x+1))
  cells<-colnames(alldata_log)
  genes<-rownames(alldata_log)
  
  require(irlba)
  xx<-alldata_log
  yy <- apply(xx, 1, function(x){
    x - mean(x)
  })
  svd.base <- irlba(yy, nv = nPC)
  centers <- apply(xx, 1, function(x){
    mean(x)
  })
  PCs <- svd.base$u%*%diag(svd.base$d[1:nPC])
  Loadings <- svd.base$v
  PCA <- list(PCs = PCs, Loadings = Loadings, Centers = centers, batch.id.forPC=batch.id, Raw=yy, cells=cells, genes=genes)
  return(PCA)
}