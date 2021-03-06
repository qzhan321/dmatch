% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projection_visualization.R
\name{projection_visualization}
\alias{projection_visualization}
\title{projection_visualization}
\usage{
projection_visualization(
  object,
  filename = NULL,
  TopCellLineNumber = 5,
  cor.threshold = 0.3,
  ShowCellNumber = 10,
  dist.method = "euclidean",
  hclust.method = "ward.D"
)
}
\arguments{
\item{object}{A dmatch class object.}

\item{filename}{The path and name of the output png file. Default is NULL, and a heatmap is printed on the screen directly.}

\item{TopCellLineNumber}{The number of primary cell lines to be kept which are highly correlated with any cell in the samples. Set the Pearson Correlation coefficients between the rest of primary cell lines with cells in the samples to zero.}

\item{cor.threshold}{If the correlation between some cells in the samples and the last TopCellLine in the reference is lower than the threshold, a warning will be generated.}

\item{ShowCellNumber}{Include only the primary cell lines which are highly correlated with more than this amount of cells in the samples.}

\item{dist.method}{The distance metric for calculating the distance between cells.}

\item{hclust.method}{The agglomeration method used by hclust() in the projection_visualization step for clustering data}
}
\value{
A dmatch class object which have slots storing raw.data, batch.id, PCA, and more information. Specifically, projection.visualization slot stores information for the weight matrix between cells in the samples and primary cell lines in the reference with enforeced sparsity, and a heat map.
}
\description{
Visualize the results of projection step and inspect the major cell types in the samples and the number of shared cell types between the samples. Sparsity is introduced in the correlation matrix by only retaining the top several primary cell lines (all others are set to zero). A minimum threshold for the correlations is provided via cor.threshold parameter. If the overall correlation between cells in the samples and the last TopCellLine in the reference panel is less than the threshold, a warning will be generated and those cells will be removed from visualization. If the overall correlations are very low, our method does not work well with those datasets.
}
\author{
Mengjie Chen, Qi Zhan
}
