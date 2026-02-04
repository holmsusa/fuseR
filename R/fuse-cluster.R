#' Perform Hierarchical Clustering on Methylation Data
#'
#' @description Produces a hierarchical clustering tree based on the input matrices of counts.
#' @param K0 Integer or numeric matrix with unmethylated counts
#' @param K1 Integer or numeric matrix with methylated counts
#' @param chr (Optional) Integer vector with chromosome info
#' @param pos (Optional) Integer vector with genomic coordinate positions
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
#' set.seed(1234)
#' K0 <- matrix(
#'   rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
#'   nrow = 100, byrow = TRUE
#' )
#' K1 <- matrix(
#'   rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
#'   nrow = 100, byrow = TRUE
#' )
#' tree <- fuse.cluster(K0, K1)
#' tree
#'
#' @export
fuse.cluster <- function(K0, K1, chr = NULL, pos = NULL) {
  # Produces a hierarchical clustering tree based on the input arrays.
  # Input: matrix K0, matrix K1
  # Output: matrix

  # Checking input
  if(!all((is.matrix(K0) || methods::is(K0, "DelayedArray")),
          (is.matrix(K1) || methods::is(K1, "DelayedArray")),
          all(dim(K0) == dim(K1))))
    stop("Wrong K0/K1 format")


  # Filling the chr.idx with chromosome information, if given
  if(length(chr) > 0) {

    # Making chromosomes integers based on the chr vector
    chr.idx <- cumsum( c( TRUE, chr[-1] != chr[-length(chr)] ) )

  } else {

    # Otherwise "one" chromosome and 1,2,3,... positions
    chr.idx <- rep(1L, nrow(K0))
    pos <- seq_len(nrow(K0))
  }

  if(!all(is.integer(chr.idx), length(chr.idx) == nrow(K0), is.integer(pos), length(pos) == nrow(K0))) stop("Wrong input for fuse_cluster_R")

  # All set, calling the fuse_cluster_R now
  tree <- .Call('fuse_cluster_R', K0, K1, chr.idx, pos)
  tree <- .Call('sort_tree_R', tree)

  return(`colnames<-`(tree, c("m1", "m2", "logl_tot", "logl_merge", "genomic_dist")))
}

