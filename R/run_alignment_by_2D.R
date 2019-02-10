#' run_alignment_by_2D
#'
#' Correct batch effects 
#'
#'
#' @author Mengjie Chen
#' @param object A dmatch class object
#' @param K The number of PCs for correcting batch effects
#' @param NCell The smallest number of cells that a selected cluster should have. Default is 100, and recommend no less than 5% of the sample size
#' @return A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information
#' @export

run_alignment_by_2D <- function(object, quantile = 0.95, K = 30, selected = NULL,
                                Nclust = NULL, steps = 20, gra_steps = 10, NCell = 100){
  
  batch.id<-object@batch.id
  Data1<-object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[1],]
  Data2<-object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[2],]
  Labels1<-object@metadata$CellType[object@metadata$batch==batch.id[1]]
  Labels2<-object@metadata$CellType[object@metadata$batch==batch.id[2]]
  
  require(MASS)
  if(is.null(selected)){
    Labels1.tab <- table(Labels1)
    Labels2.tab <- table(Labels2)
    
    Good.Labels1 <- Labels1.tab[Labels1.tab >= NCell]
    Good.Labels2 <- Labels2.tab[Labels2.tab >= NCell]
    SharedType <- intersect(names(Good.Labels1), names(Good.Labels2))
    
    CellSum <- Labels1.tab[SharedType] + Labels2.tab[SharedType]
    SharedType <- SharedType[order(CellSum)]
  } else {
    SharedType <- selected
  }
  
  if(is.null(Nclust)){
    SharedType <- SharedType
  }  else if(length(SharedType) < Nclust){
    SharedType <- SharedType
    warnings("Specified number of cell types to use is larger than the available. Use the maximum available number")
  } else {
    SharedType <- SharedType[1:Nclust]
  }
  
  if(length(SharedType) == 0){
    stop("No common cell type detected. The alignment requires at least one cell type shared by two datasets.")
  }
  
  num <- floor(K/2)
  
  all.corrected <- NULL
  for(j in 1:num){
    inter <- (2*j-1):(2*j)
    corrected_list <- run_alignment_robust_given_labels(Data1[, inter], Data2[, inter], Labels1, Labels2, SharedType,
                                                        quantile = quantile, steps = steps, gra_steps = gra_steps)
    all.corrected <- cbind(all.corrected, corrected_list$Corrected)
  }
  
  object@run_alignment_by_2D.results <- list(RefTypes = SharedType, Original = Data1, 
                                             Reference = Data2, Corrected = all.corrected,
                                             Labels1 = Labels1, Labels2 = Labels2)
  
  return(object)
  
}





run_alignment_robust <- function(object, batch.id, quantile = 0.95, selected = NULL,
                                 Nclust = NULL, steps = 20, gra_steps = 10, NCell = 100){
  
  Data1<-object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[1],1:2]
  Data2<-object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[2],1:2]
  Labels1<-object@metadata$CellType[object@metadata$batch==batch.id[1]]
  Labels2<-object@metadata$CellType[object@metadata$batch==batch.id[2]]
  #Types1 <- unique(Labels1)
  #Types2 <- unique(Labels2)
  
  #SharedType <- intersect(Types1, Types2)
  require(MASS)
  if(is.null(selected)){
    Labels1.tab <- table(Labels1)
    Labels2.tab <- table(Labels2)
    
    Good.Labels1 <- Labels1.tab[Labels1.tab >= NCell]
    Good.Labels2 <- Labels2.tab[Labels2.tab >= NCell]
    SharedType <- intersect(names(Good.Labels1), names(Good.Labels2))
    
    CellSum <- Labels1.tab[SharedType] + Labels2.tab[SharedType]
    SharedType <- SharedType[order(CellSum)]
  } else {
    SharedType <- selected
  }
  
  if(is.null(Nclust)){
    SharedType <- SharedType
  }  else if(length(SharedType) < Nclust){
    SharedType <- SharedType
    warnings("Specified number of cell types to use is larger than the available. Use the maximum available number")
  } else {
    SharedType <- SharedType[1:Nclust]
  }
  
  if(length(SharedType) == 0){
    stop("No common cell type detected. The alignment requires at least one cell type shared by two datasets.")
  }
  
  if(length(SharedType) == 1){
    X <- Data1[Labels1 == SharedType, ]
    x1 <- Data2[Labels2 == SharedType, ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    
    X <- X[subx.flag, ]
    x1 <- x1[subx1.flag, ]
    res <- run_simple_est_theta_coordinate_descent_one_cluster_cpp(X, x1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
  }
  
  if(length(SharedType) == 2){
    X <- Data1[Labels1 == SharedType[1], ]
    x1 <- Data2[Labels2 == SharedType[1], ]
    Y <- Data1[Labels1 == SharedType[2], ]
    y1 <- Data2[Labels2 == SharedType[2], ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    suby.flag <- cov.mve(Y, quantile.used = quantile*nrow(Y))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    suby1.flag <- cov.mve(y1, quantile.used = quantile*nrow(y1))$best
    
    X <- X[subx.flag, ]
    Y <- Y[suby.flag, ]
    x1 <- x1[subx1.flag, ]
    y1 <- y1[suby1.flag, ]
    
    res <- run_simple_est_theta_coordinate_descent_cpp(X, Y, x1, y1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
  }
  
  if(length(SharedType) >= 3){
    X <- Data1[Labels1 == SharedType[1], ]
    x1 <- Data2[Labels2 == SharedType[1], ]
    Y <- Data1[Labels1 == SharedType[2], ]
    y1 <- Data2[Labels2 == SharedType[2], ]
    Z <- Data1[Labels1 == SharedType[3], ]
    z1 <- Data2[Labels2 == SharedType[3], ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    suby.flag <- cov.mve(Y, quantile.used = quantile*nrow(Y))$best
    subz.flag <- cov.mve(Z, quantile.used = quantile*nrow(Z))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    suby1.flag <- cov.mve(y1, quantile.used = quantile*nrow(y1))$best
    subz1.flag <- cov.mve(z1, quantile.used = quantile*nrow(z1))$best
    
    X <- X[subx.flag, ]
    Y <- Y[suby.flag, ]
    Z <- Z[subz.flag, ]
    x1 <- x1[subx1.flag, ]
    y1 <- y1[suby1.flag, ]
    z1 <- z1[subz1.flag, ]
    
    res <- run_simple_est_theta_coordinate_descent_three_cluster_cpp(X, Y, Z, x1, y1, z1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
    
  }
  
  object@run_alignment_robust.results <- list(RefTypes = SharedType, Reference = Data2,
                                              Labels1 = Labels1, Labels2 = Labels2,
                                              Orignial = Data1, Corrected = corrected,
                                              estA = estA, estd = estd)
  
  return(object)
  
}

run_alignment_robust_given_labels <- function(Data1, Data2, Labels1, Labels2, SharedType, quantile = 0.95, steps = 20, gra_steps = 10){
  
  if(length(SharedType) == 1){
    X <- Data1[Labels1 == SharedType, ]
    x1 <- Data2[Labels2 == SharedType, ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    
    X <- X[subx.flag, ]
    x1 <- x1[subx1.flag, ]
    res <- run_simple_est_theta_coordinate_descent_one_cluster_cpp(X, x1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
  }
  
  if(length(SharedType) == 2){
    X <- Data1[Labels1 == SharedType[1], ]
    x1 <- Data2[Labels2 == SharedType[1], ]
    Y <- Data1[Labels1 == SharedType[2], ]
    y1 <- Data2[Labels2 == SharedType[2], ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    suby.flag <- cov.mve(Y, quantile.used = quantile*nrow(Y))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    suby1.flag <- cov.mve(y1, quantile.used = quantile*nrow(y1))$best
    
    X <- X[subx.flag, ]
    Y <- Y[suby.flag, ]
    x1 <- x1[subx1.flag, ]
    y1 <- y1[suby1.flag, ]
    
    res <- run_simple_est_theta_coordinate_descent_cpp(X, Y, x1, y1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
  }
  
  if(length(SharedType) >= 3){
    X <- Data1[Labels1 == SharedType[1], ]
    x1 <- Data2[Labels2 == SharedType[1], ]
    Y <- Data1[Labels1 == SharedType[2], ]
    y1 <- Data2[Labels2 == SharedType[2], ]
    Z <- Data1[Labels1 == SharedType[3], ]
    z1 <- Data2[Labels2 == SharedType[3], ]
    
    subx.flag <- cov.mve(X, quantile.used = quantile*nrow(X))$best
    suby.flag <- cov.mve(Y, quantile.used = quantile*nrow(Y))$best
    subz.flag <- cov.mve(Z, quantile.used = quantile*nrow(Z))$best
    subx1.flag <- cov.mve(x1, quantile.used = quantile*nrow(x1))$best
    suby1.flag <- cov.mve(y1, quantile.used = quantile*nrow(y1))$best
    subz1.flag <- cov.mve(z1, quantile.used = quantile*nrow(z1))$best
    
    X <- X[subx.flag, ]
    Y <- Y[suby.flag, ]
    Z <- Z[subz.flag, ]
    x1 <- x1[subx1.flag, ]
    y1 <- y1[suby1.flag, ]
    z1 <- z1[subz1.flag, ]
    
    res <- run_simple_est_theta_coordinate_descent_three_cluster_cpp(X, Y, Z, x1, y1, z1, steps, gra_steps)
    estA <- res$A
    estd <- res$d
    corrected <- t(apply(Data1%*%estA, 1, function(z){z - estd}))
    
  }
  
  res <- list(RefTypes = SharedType, Target = Data2,
              Labels1 = Labels1, Labels2 = Labels2,
              Orignial = Data1, Corrected = corrected,
              estA = estA, estd = estd)
  
  return(res)
  
}











