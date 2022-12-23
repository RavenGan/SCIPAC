#' Perform the clustering method embedded in Seurat package.
#'
#' @param sc.dat dimension reduced scRNA-seq data whose rows for cells and columns are for PCs.
#' @param res the resolution used in seurat clustering. Default is \code{res = 2.0}.
#' @return a list contains
#' \itemize{
#' \item \code{k}, the number of clusters;
#' \item \code{ct.assignment}, a data frame with one column indicating cluster assignment. Row names of \code{ct.assignment} are cell names
#' \item \code{centers}, cluster centroids. Rows are for PCs and columns are for clusters. Each cluster centroid is calculated by taking the average value of all the cells in the cluster.
#' }
#' @importFrom dplyr %>%
#' @export

seurat.ct <- function(sc.dat, res = 2.0){
  n.PC <- ncol(sc.dat)

  seurat.dat <- Seurat::CreateSeuratObject(counts = t(sc.dat))
  new.embedding <- Seurat::CreateDimReducObject(sc.dat, assay = "RNA", key = "pca_")
  seurat.dat@reductions$pca <- new.embedding
  seurat.dat <- Seurat::FindNeighbors(seurat.dat, dims = 1:n.PC)
  seurat.dat <- Seurat::FindClusters(seurat.dat, resolution = res)

  clean_batch_cluster <- as.data.frame(Seurat::Idents(seurat.dat))
  colnames(clean_batch_cluster) <- "cluster_assignment"
  k <- length(unique(Seurat::Idents(seurat.dat)))
  clean_batch_cluster$cluster_assignment <- as.numeric(clean_batch_cluster$cluster_assignment)


  centers <- sapply(c(1:k), function(ct){
    sub.bc <- rownames(clean_batch_cluster)[which(clean_batch_cluster$cluster_assignment == ct)]
    sub.dat <- t(sc.dat)[, sub.bc, drop = FALSE]
    rowSums(sub.dat)/ncol(sub.dat)
  })

  return(list("k" = k,
              "ct.assignment" = clean_batch_cluster,
              "centers" = centers))
}

