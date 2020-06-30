#' fastSVD
#'
#' Perform SVD. The input datasets can be the raw read count matrices; set the min.cells and min.genes to filter out based on any user-defined criteria; also set the LogNormalize to be TRUE to perform global-scaling normalization method “LogNormalize”. Alternaticely, the input datasets can be the already preprocessed (filtering) and nomalized ones from other packages; set min.genes and min.cells to be NULL and Lognormalize to be FALSE. 
 
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param sample.list A list of input datasets, each of which is a read count matrix with rows as genes and columns as cells.
#' @param nPC Total number of PCs to compute and store (30 by default).
#' @param min.cells Include genes expressed in at least this many cells; recommend filtering out genes expressed in less than 5 percent cells.
#' @param min.genes Include cells which at least this many genes are expressed; recommend filtering out cells with less than 200 genes expressed.
#' @return A list consists of: PCs; Loadings; Centers; batch.id for samples while performing SVD; Raw representing the input for SVD, which is a centered, combined, log, library-size normalized matrix; cells that are included and used for SVD; genes that are included and used for SVD.
#' @export

fastSVD <- function(samples.list, nPC = 30, min.cells=NULL, min.genes=NULL, LogNormalize = T){
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
  
  if (logNormalize == T) {
    alldata_norm <- apply(alldata, 2, function(x) x/sum(x))
    alldata_libsize <- apply(alldata_norm, 2, function(x) x*1000000)
    alldata_log <- apply(alldata_libsize, 2, function(x) log2(x+1))
    cells<-colnames(alldata_log)
    genes<-rownames(alldata_log) 
  } else {
    alldata_log <- alldata
    cells<-colnames(alldata_log)
    genes<-rownames(alldata_log) 
  }
  
  
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
