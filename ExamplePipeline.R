suppressPackageStartupMessages({
  library(dmatch)
})

path1<-system.file("extdata", "HC1_PBMC.barcode_gene.csv.names.gz", package = "dmatch")
path2<-system.file("extdata", "SLE_BMMC.barcode_gene.csv.names.gz", package = "dmatch")

sample1<-read.csv(gzfile(path1), row.names = 1, header = T)
sample2<-read.csv(gzfile(path2), row.names = 1, header = T)
colnames(sample1)<-paste("sample1", colnames(sample1), sep = "_")
colnames(sample2)<-paste("sample2", colnames(sample2), sep = "_")

PCA<-fastSVD(list(sample1,sample2), min.cells = 3, min.genes = 200)
PCAPlot(PCA = PCA, PCs.to.plot = c(1,2), batchs.to.plot = c(1,2), filename = "PCAPlot.png")

samples<-CreatedmatchObject(raw.data = list(sample2,sample1), batch.id = c(2,1), PCA = PCA)
path<-system.file("extdata", "cell_atlas_ref_panel", package = "dmatch")
load(path)

samples<-projection_to_reference_panel(samples, Reference = cell.line, use.genes.threshold = 0.75)
samples<-projection_visualization(samples,hclust.method = "ward.D", filename = "projection_visualization.png")

samples<-cut_groups(samples,K=5, method = "ward.D")
samples<-select_clusters(samples, quantile = 0.98)
cells.num<-samples@select.clusters$cells.num
save(cells.num, file = "select_clusters_cells_number")
shapiro.result<-samples@select.clusters$shapiro.test.pvalue
save(shapiro.result, file = "select_clusters_shapiro_result")

samples<-run_alignment_by_2D(samples, selected = c(3,4), quantile = 0.98)
plot_corrected(samples, pcs.plot = c(1:2), celltypes.plot = c(1,2,3,4,5), filename = "plot_corrected_all.png")
plot_corrected(samples, pcs.plot = c(1:2), celltypes.plot = c(3), filename = "plot_corrected_individual_cluster.png")

samples<-move_back_to_original_space(samples)
save(samples, file = "samples_after_correction_object")

