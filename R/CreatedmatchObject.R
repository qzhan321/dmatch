setClass("dmatch", representation(raw.data="data.frame", batch.id="vector", metadata="data.frame", PCA="list", project.name="character", Projection="list", Projection.visualization="list", select.clusters="list", run_alignment_by_2D.results="list", outcome.list="list"))

#suggest put the one which has more cells as the latter one in the list for raw.data

#' CreatedmatchObject
#'
#' Create a dmatch class object
#'
#'
#' @author Mengjie Chen
#' @param raw.data A list of pairwise samples for batch effects correction
#' @param batch.id Batch.id which are used to denote those two samples in the previous fastSVD step 
#' @param project The name of the project
#' @param PCA The output from fastSVD 
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA and etc information
#' @export
CreatedmatchObject<-function(raw.data, batch.id, project = "dmatchProject", PCA=PCA) 
{ 
  
  batch.id.forPC<-PCA$batch.id.forPC
  
  raw.data1<-raw.data[[1]]
  batch.id1<-batch.id[1]
  genes<-PCA$genes
  raw.data1<-raw.data1[genes,]
  cells1<-PCA$cells[batch.id.forPC==batch.id1]
  raw.data1<-raw.data1[,cells1]
  
  raw.data2<-raw.data[[2]]
  batch.id2<-batch.id[2]
  genes<-PCA$genes
  raw.data2<-raw.data2[genes,]
  cells2<-PCA$cells[batch.id.forPC==batch.id2]
  raw.data2<-raw.data2[,cells2]
  
  raw.data<-cbind(raw.data1,raw.data2)
  batch<-as.data.frame(c(rep(batch.id1,ncol(raw.data1)), rep(batch.id2,ncol(raw.data2))))
  colnames(batch)<-"batch"
  rownames(batch)<-colnames(raw.data)
  object <- new(Class = "dmatch", raw.data = raw.data, project.name = project, PCA=PCA)
  object@metadata <- batch
  object@batch.id <- batch.id
  return(object)
}

