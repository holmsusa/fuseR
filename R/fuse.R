#' @useDynLib fuseR, .registration = TRUE
NULL

#' Perform Hierarchical Clustering on Methylation Data
#'
#' @description Produces a hierarchical clustering tree based on the input matrices of counts.
#' @param K0 Integer matrix with unmethylated counts
#' @param K1 Integer matrix with methylated counts
#' @param sort Default is true, sorts the tree in optimal order
#' @return Clustering tree as a matrix
#' @examples
#' # Example: Clustering generated data
#' K0 <- matrix(sample(c(1:200), 125, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(c(1:200), 125, replace = TRUE), ncol = 5)
#' tree <- fuse.cluster(K0, K1)
#' tree
#'
#' @export
fuse.cluster <- function(K0, K1, sort = TRUE) {
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

  # All set, calling the fuse_cluster_R now
  tree <- .Call('fuse_cluster_R', K0, K1, chr.idx, pos)

  if(sort) {
    tree <- .Call('sort_tree_R', tree)

  }

  return(`colnames<-`(tree, c("m1", "m2", "logl_tot", "logl_merge", "genomic_dist")))
}


#' Sort Clustering Tree in Optimal Order
#'
#' @description Sorts the clustering tree produced by \code{\link{fuse.cluster}}.
#' @param tree Unsorted clustering tree
#' @return Sorted tree
#' @examples
#' # Example: Sorting unsorted clustering tree
#' K0 <- matrix(sample(c(1:200), 80, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(c(1:200), 80, replace = TRUE), ncol = 5)
#' tree <- fuse.cluster(K0, K1, sort = FALSE)
#' tree.sorted <- fuse.sort.tree(tree)
#' tree
#'
#' @export
fuse.sort.tree <- function(tree) {
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
#' @return A vector indicating which cluster each element in the original data frame belonged to
#' @examples
#' # Example 1: Cutting small tree in 2 segments
#' tree <- matrix(c(
#' -1, -2,  49.53106,  49.53106,  1.14473,
#' -3, -4,  78.49604,  78.49604,  1.14473,
#' -5, -6, 147.07154, 147.07154,  1.14473,
#' 1,  2,  72.98287, 201.00997,  1.14473,
#' 4,  3, 106.38879, 454.47029,  1.14473
#' ), ncol = 5, byrow = TRUE)
#'
#' segments <- fuse.cutree(tree, 2)
#' segments
#'
#' # Example 2: Cut tree based on Bayesian Information Criterion (BIC)
#' K0 <- matrix(sample(c(1:200), 500, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(c(1:200), 500, replace = TRUE), ncol = 5)
#' tree <- fuse.cluster(K0, K1)
#' tree
#'
#' k <- c(nrow(tree):1)
#' m <- nrow(K0)
#' n <- ncol(K0)
#' bic <- numeric(nrow(tree))
#'
#' # Computing BIC for every possible number of segments
#' for(i in 1:nrow(tree)) {
#'   bic[i] <- n*log(m*n)*k[i] + 2*tree[i, 4]
#' }
#' plot(bic)
#' number_of_segments <- length(bic) - (which.min(bic)-1)
#' segments <- fuse.cutree(tree, number_of_segments)
#'
#' @export
fuse.cutree <- function (tree, k) {
  # Cuts the given tree into k clusters.
  # Input: matrix tree, numeric integer k
  # Output: integer vector
  return(.Call('cuttree_R', tree, as.integer(k)))
}








