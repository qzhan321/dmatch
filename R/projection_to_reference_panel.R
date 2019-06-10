#' projection_to_reference_panel
#'
#' Project the pairwise to the reference panel to identify separate clusters in the samples
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object
#' @param Reference The reference panel
#' @param CorMethod A character string indicating which correlation coefficient (or covariance) is to be computed. Default is "pearson"
#' @param use.genes.threshold The threshold for highly variable genes used for calculating the Pearson Correlation between cells in the samples and the primary cell lines in the reference panel
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, Projection slot stores information for the correlation matrix of cells in the samples and cells in the reference panel. The correlation matrix is calculated using the common genes in the samples and the reference. Cells in the samples which have zero expression across those common genes will be filtered out. The new batch.id for the remaining cells in the reference is in batch.id.update 
#' @export
projection_to_reference_panel <- function(object, Reference, CorMethod = "pearson", use.genes.threshold=0.75){
  data<-object@raw.data
  data_norm<-apply(data,2,function(x) x/(sum(x)+0.0001))
  data_libsize<-apply(data_norm,2,function(x) x*1000000)
  
  if (!is.null(use.genes.threshold)) {
    gene.expr.mean<-apply(data_libsize,1,function(x) mean(x))
    gene.expr.var<-apply(data_libsize,1,function(x) var(x))
    gene.expr.LogVMR<-log2(gene.expr.var/(gene.expr.mean+0.0001)+1)
    gene.expr.LogVMR<-as.data.frame(gene.expr.LogVMR)
    colnames(gene.expr.LogVMR)<-c("LogVMR")
    rownames(gene.expr.LogVMR)<-rownames(data_libsize)
    data_final<-cbind(data_libsize, gene.expr.LogVMR)
    gene_LogVMR_threshold<-quantile(data_final[,ncol(data_final)], use.genes.threshold)
    data_final_filter_LogVMR<-data_final[data_final[,(ncol(data_final))]>=gene_LogVMR_threshold,]
    genes.use<-rownames(data_final_filter_LogVMR)
    Data.use<-data_libsize[genes.use,]
  } else {
    Data.use<-data_libsize
  }
  
  gene.id <- rownames(Data.use)
  gene.ref <- rownames(Reference)
  common <- intersect(gene.id, gene.ref)
  
  #check whether the Reference and Data.use subsetted by common genes produce cells whose total expression is 0. 
  #This will cause error in the following Pearson Correlation calculation step.
  a<-Data.use[common,]
  b<-Reference[common,]
  aa<-apply(a,2,function(x) sum(x))
  bb<-apply(b,2,function(x) sum(x))
  if (any(aa==0) | any(bb==0)) {
    cat("Using the intersection genes of samples and Reference cause some cells with zero total expression, remove those cells...\n")
    Data.use<-Data.use[,!(aa==0)]
    Reference<-Reference[,!(bb==0)]
    batch.id.update<-object@metadata$batch[!(aa==0)]
  }
  
  logxx <- apply(Data.use[common, ], 2, function(x){
    log(x+0.1)
  })
  
  selected.cell.line <- apply(Reference[common, ], 2, function(x){
    x - mean(x)
  })
  
  n <- ncol(logxx)
  m <- ncol(Reference)
  
  cor.mat <- matrix(0, nrow = n, ncol = m)
  
  for(j in 1:m){
    for(i in 1:n){
      cor.mat[i, j] <- cor(logxx[, i], selected.cell.line[, j], method = CorMethod)
    }
  }
  rownames(cor.mat) <- colnames(Data.use)
  colnames(cor.mat) <- colnames(Reference)
  ReferenceNames <- colnames(Reference) 
  SampleNames <- colnames(Data.use)
  object@Projection <- list("cor.mat"=cor.mat, "ReferenceNames"=ReferenceNames, "SampleNames"=SampleNames, "batch.id.update"=batch.id.update)
  return(object)
}

