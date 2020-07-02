setClass("dmatch", representation(pairwiseSamples="data.frame", batch.id="vector", metadata="data.frame", PCA="list", project.name="character", Projection="list", Projection.visualization="list", cut_groups="list", select.clusters="list", batch_effects_vector_angles="matrix", run_alignment_by_2D.results="list", outcome.list="list"))


#' CreatedmatchObject
#'
#' Create a dmatch class object. It is common to correct batch effects for more than two samples. We can run fastSVD with all those samples (as a list) at once. From this step on, we correct batch effects for pairwise samples at a time. We recommend to use the same dataset (sample) which has the most cells/the best quality as the reference sample, and use it to align the rest samples one at a time. The reference sample will be the latter in the pairwiseSamples.list.

#' @author Mengjie Chen, Qi Zhan
#' @param pairwiseSamples.list Same as the input pairwiseSamples.list, or subset of the pairwiseSamples.list, for fastSVD. A list of pairwise samples for batch effects correction, each of which is a raw read count dataset. Alternatively, those samples can be the preprocessed and Lognormalized ones; then need to specify no preprocessing and no LogNormalize in the next step projection_to_reference_panel.
#' @param batch.id Batch.ids which are used to denote those two samples in the previous fastSVD step.
#' @param project Project name (string).
#' @param PCA The output from fastSVD, i.e., the list returned by fastSVD function.
#' @return Initializes the dmatch object. Return a dmatch class object which have slots storing pairwiseSamples.list, batch.id, PCA, and etc. 
#' @export

CreatedmatchObject<-function(pairwiseSamples.list, batch.id, project = "dmatchProject", PCA=PCA) 
{ 
  sample1<-pairwiseSamples.list[[1]]
  sample1.batch<-batch.id[1]
  
  sample2<-pairwiseSamples.list[[2]]
  sample2.batch<-batch.id[2]
  
  samples<-cbind(sample1,sample2)
  batch<-as.data.frame(c(rep(sample1.batch,ncol(sample1)), rep(sample2.batch,ncol(sample2))))
  colnames(batch)<-"batch"
  rownames(batch)<-colnames(samples)
  object <- new(Class = "dmatch", pairwiseSamples = samples, project.name = project, PCA=PCA)
  object@metadata <- batch
  object@batch.id <- batch.id
  
  return(object)
}

