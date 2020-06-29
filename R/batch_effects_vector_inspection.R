#' batch_effects_vector_inspection
#'
#' Function for calculating batch effects vectors and their angles, which can aid determining the right amount of PCs to correct for batch effects.
#' Batch effects vectors should be understood in the following way. When we correct batch effects for a pair of samples in 2-D PC space (PC1 and 2, PC3 and 4, ...), assume there are four shared clusters (1,2,3,4) in each of those two samples. We can draw a line (with direction) connecting the center of cluster 1 in the reference sample and the center of cluster 1 in the to-be-corrected sample. Direction of the line can be understood in the following way: positive means that the center of cluster 1 in the to-be-corrected reference sample is to the right of the center of cluster 1 in the reference sample; negative means the other way. Such lines are batch effects vectors. 
#' This function calculates the angles between those batch effects vectors. If the angles are small, for example, if the angle between batch effect vector for cluster 1 and that for cluster 2 is small, it means that the adjustments/correction for cluster 1 and cluster 2 are fairly similar. On the other hand, if the angle between batch effect vector for cluster 1 and that for cluster 2 is big, for example, almost 180 degrees, then whatever the adjustment/correction required for cluster 1 will be opposite to whatever adjustment/correction required for cluster 2. So the corrections required for different clusters are quite different.
#' Often, the default number of PCs (top 30 PCs) works well for batch effects correction. However, we can use this function to make sure that whether the angles between different batch effects vectors are overall acceptable (small or moderate sizes). Sometimes the batch effects are consistent only for top 6 PCs for example. Then we will correct batch effects on top 6 PCs.   
#' 
#'
#' @author Mengjie Chen, Qi Zhan
#' @param object A dmatch class object
#' @param quantile The minimum number of the data points regarded as good points
#' @param PCs_to_check Which PCs to check, for example, if set to be 1:10, this function will check the batch effects angles in every 2D PC space (PC1 and 2, PC3 and 4, ..., PC9 and 10)
#' @param clusters_to_check Check the angles between which clusters, more specifically, the candidate clusters for anchors. For example, if set to be list(c(1,2),c(2,3)), this function will check the angles between batch effects for cluster 1 and 2, and those for cluster 2 and 3, in the specified PC space
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, batch_effects_vector_angles slot stores information for the angles between batch effects vectors for different anchor cluster
#' @export

batch_effects_vector_inspection <- function(object, quantile=0.95, PCs_to_check, clusters_to_check){
  
  if (is.null(PCs_to_check)) {
    stop("`PCs_to_check` is not specified")
  }
  
  if (is.null(clusters_to_check)) {
    stop("`clusters_to_check` is not specified")
  }
  
  PCA<-object@PCA
  batch.id<-object@batch.id
  batch.id.forPC<-PCA$batch.id.forPC
  batch.id.update<-object@Projection$batch.id.update
  batch.id.update1<-batch.id.update[batch.id.update==batch.id[1]]
  batch.id.update2<-batch.id.update[batch.id.update==batch.id[2]]
  
  Data1 <- PCA$PCs[names(batch.id.update1),]
  Data2 <- PCA$PCs[names(batch.id.update2),]
  
  Labels1 <- object@cut_groups$CellType[object@Projection$batch.id.update==batch.id[1]]
  Labels2 <- object@cut_groups$CellType[object@Projection$batch.id.update==batch.id[2]]
  
  require(MASS)
  batch_effects_vector_angles<-matrix(0,nrow = num,ncol = length(clusters_to_check))
  num <- floor(PCs_to_check/2)
  all.corrected <- NULL
  for (j in 1:num) {
    inter <- (2*j-1):(2*j)
    for (k in 1:length(clusters_to_check)) {
      cluster1<-clusters_to_check[[k]][1]
      cluster2<-clusters_to_check[[k]][2]
      X <- Data1[Labels1 == cluster1, inter]
      x1 <- Data2[Labels2 == cluster1, inter]
      subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
      subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
      X <- X[subx.flag, ]
      x1 <- x1[subx1.flag, ]
      res <- run_simple_no_rotation_one_cluster_cpp(X, x1)
      #estA <- res$A
      estd1 <- res$d
      
      X <- Data1[Labels1 == cluster2, inter]
      x1 <- Data2[Labels2 == cluster2, inter]
      subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
      subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
      X <- X[subx.flag, ]
      x1 <- x1[subx1.flag, ]
      res <- run_simple_no_rotation_one_cluster_cpp(X, x1)
      #estA <- res$A
      estd2 <- res$d
      
      cluster1.vector<-estd1
      cluster2.vector<-estd2
      cluster12_theta <- acos( sum(cluster1.vector*cluster2.vector) / ( sqrt(sum(cluster1.vector*cluster1.vector)) * sqrt(sum(cluster2.vector*cluster2.vector)) ) )
      cluster12_theta <- cluster12_theta*180/pi
      
      batch_effects_vector_angles[j,k] <- cluster12_theta
      
      colnames(batch_effects_vector_angles)[k] <- paste0("clusters", cluster1, cluster2) 
    }
    rownames(batch_effects_vector_angles)[j] <- paste0("PC",inter[1],inter[2])
  }
  object@batch_effects_vector_angles <- batch_effects_vector_angles
    
}