% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.sc.bulk.dat.R
\name{obtain.preprocessed.data}
\alias{obtain.preprocessed.data}
\title{Normalize scRNA-seq data and find highly variable genes without scaling using Seurat}
\usage{
obtain.preprocessed.data(exprs.data, hvg = 1000)
}
\arguments{
\item{exprs.data}{single-cell data matrix. Rows stand for genes and columns stand for cells.}

\item{hvg}{the number of highly variable genes. The default is 1000.}
}
\value{
a single-cell data matrix. Rows stand for genes and columns stand for cells.
}
\description{
Normalize scRNA-seq data and find highly variable genes without scaling using Seurat
}
