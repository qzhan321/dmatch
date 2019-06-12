#' select_clusters
#'
#' Calculate the shapiro test for normality for each cluster in each sample, and return the table for the cell numbers for each cluster in each sample
#'
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object
#' @param quantile The minimum number of the data points regarded as good points
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, select.clusters slot stores information for shapiro test result and numbers of cells in each cluster in each sample
#' @export
select_clusters<-function(object, quantile=0.95) {
  require(MASS)
  PCA<-object@PCA
  batch.id<-object@batch.id
  batch.id.forPC<-PCA$batch.id.forPC
  batch.id.update<-object@Projection$batch.id.update
  batch.id.update1<-batch.id.update[batch.id.update==batch.id[1]]
  batch.id.update2<-batch.id.update[batch.id.update==batch.id[2]]
  
  Data1 <- PCA$PCs[names(batch.id.update1),]
  Data2 <- PCA$PCs[names(batch.id.update2),]
  
  Labels1 <- object@Projection$CellType[object@Projection$batch.id.update==batch.id[1]]
  Labels2 <- object@Projection$CellType[object@Projection$batch.id.update==batch.id[2]]
  
  #Shapiro-Wilk test
  shapiros1<-NULL
  for (k in 1:10) {
    Data1_PC<-Data1[,k]
    for (j in 1:length(unique(Labels1))) {
      if (length(Data1_PC[Labels1==j])<3) {
        shapiro1<-list("p.value"=NA)
      } else {
        mm<-Data1_PC[Labels1==j]
        if (quantile*nrow(as.data.frame(mm)) >= 17) {
          subx.flag <- cov.mve(as.data.frame(mm), quantile.used = round(quantile*nrow(as.data.frame(mm))))$best
          nn<-mm[subx.flag]
          nn<-(nn-mean(nn))/sd(nn)
          shapiro1<-shapiro.test(nn) 
        } else {
          mm<-(nn-mean(mm))/sd(mm)
          shapiro1<-shapiro.test(mm) 
        }
      }
      shapiros1<-c(shapiros1,shapiro1$p.value)
    }
    
  }
  nrow<-length(unique(Labels1))
  pvalue<-matrix(shapiros1, nrow = nrow, ncol=10, byrow = F)
  rownames(pvalue) <- paste0("cluster", seq_len(nrow))
  colnames(pvalue) <- paste0("PC", seq_len(10))
  a<-as.data.frame(table(Labels1))
  colnames(a)<-c("Labels","number_of_cells")
  rownames(a)<-rownames(pvalue)
  pvalue<-cbind(a,pvalue)
  
  mean_shapiro1<-apply(pvalue[,3:12], 1, function(x) sum(-log(x))/10)
  
  shapiros2<-NULL
  for (k in 1:10) {
    Data2_PC<-Data2[,k]
    for (j in 1:length(unique(Labels2))) {
      if (length(Data2_PC[Labels2==j])<3) {
        shapiro2<-list("p.value"=NA)
      } else {
        mm<-Data2_PC[Labels2==j]
        if (quantile*nrow(as.data.frame(mm)) >= 17) {
          subx.flag <- cov.mve(as.data.frame(mm), quantile.used = round(quantile*nrow(as.data.frame(mm))))$best
          nn<-mm[subx.flag]
          nn<-(nn-mean(nn))/sd(nn)
          shapiro2<-shapiro.test(nn) 
        } else {
          mm<-(mm-mean(mm))/sd(mm)
          shapiro2<-shapiro.test(mm) 
        }
      }
      shapiros2<-c(shapiros2,shapiro2$p.value)
    }
    
  }  
  nrow<-length(unique(Labels2))
  pvalue<-matrix(shapiros2, nrow = nrow, ncol=10, byrow = F)
  rownames(pvalue) <- paste0("cluster", seq_len(nrow))
  colnames(pvalue) <- paste0("PC", seq_len(10))
  a<-as.data.frame(table(Labels2))
  colnames(a)<-c("Labels","number_of_cells")
  rownames(a)<-rownames(pvalue)
  pvalue<-cbind(a,pvalue)
  
  mean_shapiro2<-apply(pvalue[,3:12], 1, function(x) sum(-log(x))/10)
  
  shapiro.test.pvalue<-matrix(c(mean_shapiro1,mean_shapiro2), nrow = length(mean_shapiro1), ncol = 2, byrow = F)
  rownames(shapiro.test.pvalue)<-paste("cell cluster", seq_along(1:nrow(shapiro.test.pvalue)), sep="_")
  colnames(shapiro.test.pvalue)<-paste("batch", batch.id[seq_along(batch.id)], sep = "_")
  
  xx<-as.data.frame(table(Labels1))
  xx<-xx[,-1, drop=F]
  yy<-as.data.frame(table(Labels2))
  yy<-yy[,-1, drop=F]
  cells.num<-cbind(xx,yy)
  rownames(cells.num)<-paste("cell cluster", seq_along(1:nrow(cells.num)), sep="_")
  colnames(cells.num)<-paste("batch", batch.id[seq_along(batch.id)], sep = "_")
  
  object@select.clusters<-list("shapiro.test.pvalue"=shapiro.test.pvalue, "cells.num"=cells.num)
  
  return(object)  
}
