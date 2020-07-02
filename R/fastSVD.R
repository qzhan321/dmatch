#' fastSVD
#'
#' Perform SVD. The input datasets are the already preprocessed (filtering) and global-scale-nomalized datasets from Seurat. 
 
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param samples.list A list of input datasets preprocessed by Seurat.
#' @param nPC Total number of PCs to compute and store (30 by default).
#' @return A list consists of: PCs; Loadings; Centers (centers of samples before SVD); batch.id.forPC (batch id of cells in the samples for SVD); Raw representing the input for SVD; cells are the cell names for cells that are included and used for SVD; genes are the gene names for genes that are included and used for SVD.
#' @export

fastSVD <- function(samples.list, nPC = 30){
  
  samples <- samples.list[[1]]
  batch.id <- rep(1,ncol(samples))
  for (i in 2:length(samples.list)) {
    sample <- samples.list[[i]]
    samples <- cbind(samples, sample)
    batch.id <- c(batch.id, rep(i, ncol(sample)))
  }  
  
  names(batch.id)<-colnames(samples)
  cells<-colnames(samples)
  genes<-rownames(samples) 
  
  require(irlba)
  # center samples before SVD
  xx<-samples
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
