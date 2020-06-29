setClass("dmatch", representation(pairwiseSample.list="data.frame", batch.id="vector", metadata="data.frame", PCA="list", project.name="character", Projection="list", Projection.visualization="list", cut_groups="list", select.clusters="list", batch_effects_vector_angles="matrix", run_alignment_by_2D.results="list", outcome.list="list"))


#' CreatedmatchObject
#'
#' Create a dmatch class object. It is common to correct batch effects for more than two samples. We can run fastSVD with all those samples (as a list) at once. From this step on, we correct batch effects for pairwise samples at a time. We recommend to use the same dataset (sample) which has the most cells as the reference sample, and use this same reference sample to align the rest samples one at a time. The reference sample will be the latter in the pairwisepairwiseSample.list.

#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param pairwiseSample.list Same as the input pairwiseSample.list, or subset of the pairwiseSample.list, for fastSVD. A list of pairwise samples for batch effects correction, each of which is a raw read count dataset. Alternatively, those samples can be the preprocessed and Lognormalized ones; then need to specify no preprocessing and no LogNormalize in the next step projection_to_reference_panel.
#' @param batch.id Batch.ids which are used to denote those two samples in the previous fastSVD step.
#' @param project Project name (string).
#' @param PCA The output from fastSVD, i.e., the list returned by fastSVD function.
#' @return Initializes the dmatch object. Return a dmatch class object which have slots storing pairwiseSample.list, batch.id, PCA, and etc. This step only includes the cells and genes used in fastSVD step for the samples provided in the pairwiseSample.list.
#' @export
CreatedmatchObject<-function(pairwiseSample.list, batch.id, project = "dmatchProject", PCA=PCA) 
{ 
  
  batch.id.forPC<-PCA$batch.id.forPC
  
  sample1<-pairwiseSample.list[[1]]
  sample1.batch<-batch.id[1]
  genes<-PCA$genes
  sample1<-sample1[genes,]
  sample1.cells<-PCA$cells[batch.id.forPC==sample1.batch]
  sample1<-sample1[,sample1.cells]
  
  sample2<-pairwiseSample.list[[2]]
  sample2.batch<-batch.id[2]
  genes<-PCA$genes
  sample2<-sample2[genes,]
  sample2.cells<-PCA$cells[batch.id.forPC==sample2.batch]
  sample2<-sample2[,sample2.cells]
  
  pairwiseSample.list.final<-cbind(sample1,sample2)
  batch<-as.data.frame(c(rep(sample1.batch,ncol(sample1)), rep(sample2.batch,ncol(sample2))))
  colnames(batch)<-"batch"
  rownames(batch)<-colnames(pairwiseSample.list.final)
  object <- new(Class = "dmatch", pairwiseSample.list = pairwiseSample.list.final, project.name = project, PCA=PCA)
  object@metadata <- batch
  object@batch.id <- batch.id
  
  return(object)
}

