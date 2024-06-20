
#' Normalize scRNA-seq data and find highly variable genes without scaling using Seurat
#'
#' @param exprs.data single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @param hvg the number of highly variable genes. The default is 1000.
#'
#' @return a single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @importFrom dplyr %>%
obtain.preprocessed.data <- function(exprs.data, hvg = 1000){
  exprs.data <- Seurat::CreateSeuratObject(exprs.data, project = "CreateSeuratObject", min.cells = 3, min.features = 200)
  exprs.data <- Seurat::NormalizeData(exprs.data, normalization.method = "LogNormalize", scale.factor = 10000)
  exprs.data <- Seurat::FindVariableFeatures(exprs.data, selection.method = "vst", nfeatures = hvg)

  exprs.data.norm <- exprs.data[["RNA"]]$data
  exprs.data.variable.genes <- Seurat::VariableFeatures(exprs.data)

  processed.data <- exprs.data.norm[exprs.data.variable.genes, ] %>%
    as.matrix()

  return(processed.data)
}



#' Pre-process both single-cell and bulk data
#'
#' This function first finds overlap genes between scRNA-seq and bulk data,
#' then normalizes scRNA-seq data and find variable genes,
#' log-transforms bulk RNA-seq data,
#' and uses high variable genes for both scRNA-seq and bulk data
#'
#' @param sc.dat single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @param bulk.dat bulk data matrix. Rows stand for genes and columns stand for cells.
#' @param hvg the number of highly variable genes. The default is \code{1000}.
#'
#' @return A list object contains pre-processed single-cell and bulk data
#' @importFrom dplyr %>%
#' @export

preprocess.sc.bulk.dat <- function(sc.dat, bulk.dat, hvg = 1000){
  overlap.genes <- intersect(rownames(sc.dat), rownames(bulk.dat))

  sc.dat.new <- sc.dat[overlap.genes, ]

  bulk.dat.new <- bulk.dat[overlap.genes, ] %>% as.matrix()
  bulk.dat.new <- log(bulk.dat.new + 1)


  sc.dat.preprocessed <- obtain.preprocessed.data(exprs.data = sc.dat.new, hvg = hvg)
  bulk.dat.preprocessed <- bulk.dat.new[rownames(sc.dat.preprocessed), ]

  return(list("sc.dat.preprocessed" = sc.dat.preprocessed,
              "bulk.dat.preprocessed" = bulk.dat.preprocessed))
}
