
fastSVD <- function(samples.list, nPC = 30, min.cells=NULL, min.genes=NULL){
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
  
  alldata_norm <- apply(alldata, 2, function(x) x/sum(x))
  alldata_libsize <- apply(alldata_norm, 2, function(x) x*1000000)
  alldata_log <- apply(alldata_libsize, 2, function(x) log(x+1))
  cells<-colnames(alldata_log)
  genes<-rownames(alldata_log)
  
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
  PCA <- list(PCs = PCs, Loadings = Loadings, Centers = centers, batch.id.forPC=batch.id, Raw=yy, cells=cells, genes=genes)
  return(PCA)
}



PCAPlot<-function(PCA, PCs.to.plot, filename) {
  i<-PCs.to.plot[1]
  j<-PCs.to.plot[2]
  PC1 <- PCA$PCs[, i]
  PC2 <- PCA$PCs[, j]
  batch.id.forPC<-PCA$batch.id.forPC
  png(filename, res = 400, height = 8, width = 8, unit = "in")
  plot(PC1, PC2, type = "n")
  text(PC1, PC2, batch.id.forPC, col = batch.id.forPC, cex = 0.8)
  dev.off()
}






setClass("dmatch", representation(raw.data="data.frame", batch.id="vector", metadata="data.frame", PCA="list", project.name="character", Projection="list", Projection.visualization="list", select.clusters="list", run_alignment_by_2D.results="list", outcome.list="list"))

#suggest put the one which has more cells as the latter one in the list for raw.data


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





#two references, decide which one to use

projection_to_reference_panel <- function(object, Reference, CorMethod = "pearson", use.genes.threshold=NULL){
  data<-object@raw.data
  data_norm<-apply(data,2,function(x) x/sum(x))
  data_libsize<-apply(data_norm,2,function(x) x*1000000)
  
  if (!is.null(use.genes.threshold)) {
    gene.expr.mean<-apply(data_libsize,1,function(x) mean(x))
    gene.expr.var<-apply(data_libsize,1,function(x) var(x))
    gene.expr.LogVMR<-log(gene.expr.var/(gene.expr.mean+0.0001)+1)
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
  object@Projection <- list("cor.mat"=cor.mat, "ReferenceNames"=ReferenceNames, "SampleNames"=SampleNames)
  return(object)
}


projection_visualization <- function(object, filename = NULL, TopCellLineNumber = 5, ShowCellNumber = 20, dist.method="euclidean", hclust.method="ward.D"){
  
  Projected <- object@Projection$cor.mat
  ReferenceNames <- object@Projection$ReferenceNames
  SampleNames <- object@Projection$SampleNames
  
  weights.mat <- apply(Projected, 1, function(x){
    y <- x
    cutoff <- x[order(x, decreasing=T)[TopCellLineNumber]]
    y[which(x<cutoff)] <- 0
    y
  })
  
  flag <- apply(weights.mat, 1, function(x){
    length(x[x!=0]) >= ShowCellNumber
  })
  
  kkk <- as.matrix(weights.mat)
  
  require(gplots)
  palette.gr.marray2 <- colorRampPalette(c("white", "red"))(16)
  colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki",
                 "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",
                 "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
  
  png(filename, res = 400, height = 8, width = 8, unit = "in")
  dist.method<-dist.method
  hclust.method<-hclust.method
  aa <- heatmap.2(kkk[flag, ], trace = "none", col = palette.gr.marray2, symbreaks = F,
                  labRow = ReferenceNames[flag], labCol = NA,  ColSideColors = colorlist[as.numeric(object@metadata$batch)],
                  key = F, margins = c(8, 15), distfun=function(x) dist(x,method = dist.method), hclustfun=function(x) hclust(x,method= hclust.method))
  dev.off()
  # bb <- rev(aa$colInd)
  
  bbb = kkk[flag, ]
  rownames(bbb) = ReferenceNames[flag]
  colnames(bbb) = SampleNames
  object@Projection.visualization<-list("Weight.max" = bbb, "ProjectionHeatmap" = aa)
 
  return(object)
  
}


cut_groups <- function(object, K, method="ward.D"){
  
  WeightMat <- object@Projection.visualization$Weight.max
  dd <- dist(t(WeightMat))
  ee <- hclust(dd, method)
  groups <- cutree(ee, k=K)
  groups<-as.data.frame(groups)
  rownames(groups)<-colnames(object@raw.data)
  colnames(groups)<-"CellType"
  object@metadata<-cbind(object@metadata,groups) 
  return(object)
}




select_clusters<-function(object, quantile=0.95) {
  require(MASS)
  PCA<-object@PCA
  batch.id<-object@batch.id
  batch.id.forPC<-PCA$batch.id.forPC
  Data1 <- PCA$PCs[batch.id.forPC==batch.id[1],]
  Data2 <- PCA$PCs[batch.id.forPC==batch.id[2],]
  Labels1 <- object@metadata$CellType[object@metadata$batch==batch.id[1]]
  Labels2 <- object@metadata$CellType[object@metadata$batch==batch.id[2]]
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


plot_corrected <- function(object, filename = NULL, pcs.plot, celltypes.plot){
  
  SharedType <- object@run_alignment_by_2D.results$RefTypes
  
  Data2 <- object@run_alignment_by_2D.results$Reference
  Data1 <- object@run_alignment_by_2D.results$Original
  
  Labels1 <- object@run_alignment_by_2D.results$Labels1
  Labels2 <- object@run_alignment_by_2D.results$Labels2
  
  Corrected <- object@run_alignment_by_2D.results$Corrected
  
  pc.num1<-pcs.plot[1]
  pc.num2<-pcs.plot[2]
  
  png(filename, res = 400, height = 6, width = 10, unit = "in")
  par(mfrow = c(1, 2))
  xlim <- c(min(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))-1, max(c(Data1[, pc.num1], Data2[, pc.num1], Corrected[, pc.num1]))+1)
  ylim <- c(min(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))-1, max(c(Data1[, pc.num2], Data2[, pc.num2], Corrected[, pc.num2]))+1)
  
  plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "Before Correction")
  
  for(i in 1:length(celltypes.plot)){
    X <- Data1[Labels1 == celltypes.plot[i], pcs.plot]
    x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
    points(x1, col = "black", pch = i)
    points(X, col = "red", pch = i)
  }
  
  plot(0, 0, xlim = xlim, ylim = ylim, col = "white", main = "After Correction")
  
  for(i in 1:length(celltypes.plot)){
    x1 <- Data2[Labels2 == celltypes.plot[i], pcs.plot]
    X <- Corrected[Labels1 == celltypes.plot[i], pcs.plot]
    points(x1, col = "black", pch = i)
    points(X, col = "green", pch = i)
  }
  dev.off()
  
}





move_back_to_original_space <- function(object){
  batch.id<-object@batch.id
  Final1 <- t(object@PCA$Raw[object@PCA$batch.id.forPC==batch.id[1],]) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Corrected - object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[1],])
  Final2 <- t(object@PCA$Raw[object@PCA$batch.id.forPC==batch.id[2],]) + object@PCA$Loadings%*%t(object@run_alignment_by_2D.results$Corrected - object@PCA$PCs[object@PCA$batch.id.forPC==batch.id[2],])
  object@outcome.list<-list(Final1,Final2)
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











