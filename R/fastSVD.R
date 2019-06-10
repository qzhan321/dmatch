#' fastSVD
#'
#' Perform SVD
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param sample.list A list of raw inputs, each of which is a read count matrix (non-normalized) with rows as genes and columns as cells
#' @param nPC Total number of PCs to compute and store (30 by default)
#' @param min.cells Include genes with detected expression in at least this many cells
#' @param min.genes Include cells where at least this many genes are detected
#' @return A list consists of: PCs; Loadings; Centers; batch.id for samples while performing SVD; Raw representing the input for SVD, which is a centered, combined, log, library-size normalized matrix; cells that are included and used for SVD; genes that are included and used for SVD
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
  
  names(batch.id)<-colnames(alldata)
  
  alldata_norm <- apply(alldata, 2, function(x) x/sum(x))
  alldata_libsize <- apply(alldata_norm, 2, function(x) x*1000000)
  alldata_log <- apply(alldata_libsize, 2, function(x) log2(x+1))
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
  rownames(PCs)<-cells
  colnames(PCs)<-paste0("PC", seq_along(1:nPC))
  rownames(Loadings)<-genes
  colnames(Loadings)<-paste0("PC", seq_along(1:nPC))
  PCA <- list(PCs = PCs, Loadings = Loadings, Centers = centers, batch.id.forPC=batch.id, Raw=yy, cells=cells, genes=genes)
  return(PCA)
}
