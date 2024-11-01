#dyn.unload('./src/test/C/island.so')
#dyn.load('./src/C/island.so')

#' @useDynLib island, .registration = TRUE
NULL

#' Perform Hierarchical Clustering on Methylation Data
#'
#' @description Produces a hierarchical clustering tree based on the input matrices of counts.
#' @param K0 Integer matrix with unmethylated counts
#' @param K1 Integer matrix with methylated counts
#' @return Unsorted clustering tree as a matrix
#' @export
island.cluster <- function(K0, K1) {
  # Produces a hierarchical clustering tree based on the input arrays.
  # Input: matrix K0, matrix K1
  # Output: matrix 
  
  # Checking input
  stopifnot(is.matrix(K0), is.matrix(K1), dim(K0) == dim(K1))

  # Initializing chromosome vector (chr.idx) and genomic location vector (pos)
  chr.idx <- as.integer(rep(0, nrow(K0)))
  pos <- as.integer(rep(0, nrow(K0)))
  
  # Saving the chromosomes into vector 
  tmp <- rownames(K0)
  chr <- sub( '^(.*)[.](\\d+)$', '\\1', tmp )
  
  # Filling the chr.idx with chromosome information, if there is any
  if(length(chr) > 0) {

    # Making chromosomes integers based on the chr vector
    chr.idx <- cumsum( c( TRUE, chr[-1] != chr[-length(chr)] ) )

    # Extracting the genomic locations into the pos vector
    pos <- as.integer( sub( '^(.*)[.](\\d+)$', '\\2', tmp ) )
  }
  
  stopifnot(is.integer(chr.idx), length(chr.idx) == nrow(K0), is.integer(pos), length(pos) == nrow(K0))

  # All set, calling the island_cluster_R now
  return (.Call('island_cluster_R', K0, K1, chr.idx, pos))
}
	

#' Sort Clustering Tree in Optimal Order
#'
#' @description Sorts the clustering tree produced by \code{\link{island.cluster}}.
#' @param tree Unsorted clustering tree
#' @return Sorted tree
#' @export
island.sort.tree <- function(tree) {
  # Sorts the tree in optimal order, returns 
  # Input: matrix tree
  # Output: matrix 
  return (.Call('sort_tree_R', tree))
}

#' Cut Hierarchical Clustering Tree into Clusters
#'
#' @description Divides the clustering tree into a specified number of clusters.
#' @param tree Sorted clustering tree
#' @param k Number of clusters
#' @return A vector indicating which cluster each element in the original dataframe belonged to 
#' @export
island.cutree <- function (tree, k) {
  # Cuts the given tree into k clusters. 
  # Input: matrix tree, numeric integer k
  # Output: integer vector 
  return(.Call('cuttree_R', tree, as.integer(k)))
}








