#' @useDynLib fuseR, .registration = TRUE
NULL

#' Perform Hierarchical Clustering on Methylation Data
#'
#' @description Produces a hierarchical clustering tree based on the input matrices of counts.
#' @param K0 Integer or numeric matrix with unmethylated counts
#' @param K1 Integer or numeric matrix with methylated counts
#' @param sort Default is true, sorts the tree in optimal order
#'
#' @return
#' Clustering tree as a matrix, with the same structure as hclust, including the columns
#' \describe{
#'   \item{m1}{ID of first data point. <0 is original point, >0 is row of previous merge.}
#'   \item{m2}{ID of second data point. <0 is original point, >0 is row of previous merge.}
#'   \item{logl_tot}{Change in total log-likelihood of merge.}
#'   \item{logl_merge}{Total cost of points in this merge.}
#'   \item{genomic_dist}{Penalty for genomic distance between points.}
#' }
#'
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
