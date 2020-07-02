#' projection_to_reference_panel
#'
#' Project the pairwise samples to the reference panel to identify separate clusters in them.
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object.
#' @param Reference The reference panel.
#' @param CorMethod A character string indicating which correlation coefficient (or covariance) is to be computed. Default is "pearson".
#' @return A dmatch class object which have slots storing pairwiseSample.list, batch.id, PCA, and more information. Specifically, Projection slot stores information for the correlation matrix of cells in the samples and cells in the reference panel. The correlation matrix is calculated using the common genes in the samples and the reference. Cells in the samples which have zero expression across those common genes will be filtered out. The new batch.id for the remaining cells in the reference is in batch.id.update. 
#' @export

projection_to_reference_panel <- function(object, Reference, CorMethod = "pearson"){
  samples<-object@pairwiseSamples
  
  gene.id <- rownames(samples)
  gene.ref <- rownames(Reference)
  common <- intersect(gene.id, gene.ref)
  cat(paste("There are in total", length(common), "common genes for the data and the reference; using those common genes to calculate the correlations between cells in the data and reference..."))
  
  #check whether the Reference and samples subsetted by common genes produce cells whose total expression is 0. 
  #This will cause error in the following Pearson Correlation calculation step.
  samples <- samples[common,]
  Reference <- Reference[common,]
  aa <- apply(samples,2,function(x) sum(x))
  bb <- apply(Reference,2,function(x) sum(x))
  if (any(aa==0) | any(bb==0)) {
    cat("Using the intersection genes of samples and Reference cause some cells with zero total expression, remove those cells...\n")
    samples<-samples[,!(aa==0)]
    Reference<-Reference[,!(bb==0)]
    batch.id.update<-object@metadata$batch[!(aa==0)]
    names(batch.id.update)<-(rownames(object@metadata))[!(aa==0)]
  } else {
    batch.id.update<-object@metadata$batch
    names(batch.id.update)<-rownames(object@metadata)
  }
  
  Reference_centered <- apply(Reference, 2, function(x){
    x - mean(x)
  })
  
  n <- ncol(samples)
  m <- ncol(Reference)
  
  cor.mat <- matrix(0, nrow = n, ncol = m)
  cor.mat <- cor(samples, Reference_centered, method = CorMethod)
  
  ReferenceNames <- colnames(Reference) 
  SampleNames <- colnames(samples)
  rownames(cor.mat) <- SampleNames 
  colnames(cor.mat) <- ReferenceNames

  object@Projection <- list("cor.mat"=cor.mat, "ReferenceNames"=ReferenceNames, "SampleNames"=SampleNames, "batch.id.Proj.To.Ref"=batch.id.update)
  return(object)
}

