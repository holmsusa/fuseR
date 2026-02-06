#' Perform Hierarchical Clustering on Methylation Data
#'
#' @description Produces a hierarchical clustering tree based on the input matrices of counts.
#' @param x Input object. One of:
#' \describe{
#'   \item{matrix}{Unmethylated count matrix (K0).}
#'   \item{BSseq}{A \code{BSseq} object.}
#'   \item{methrix}{A \code{methrix} object.}
#' }
#' @param ... K1 Integer or numeric matrix with methylated counts
#' #' @param ... Additional arguments depending on input type:
#' \describe{
#'   \item{K1}{Methylated count matrix (required if \code{x} is a matrix).}
#'   \item{chr}{(Optinal) Chromosome labels, one per CpG (matrix input only).}
#'   \item{pos}{(Optional) positions, one per CpG (matrix input only).}
#' }
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
fuse.cluster <- function(x, ...) {
  dots <- list(...)

  # If matrix input and first unnamed arg exists → treat as K1
  if (is.matrix(x) && length(dots) > 0 && is.null(names(dots)[1])) {
    dots$K1 <- dots[[1]]
    dots[[1]] <- NULL
  }

  UseMethod("fuse.cluster", x)
}


#' @rdname fuse.cluster
#' @param K1 Methylated count matrix
#' @param chr Optional chromosome vector
#' @param pos Optional position vector
#' @export
fuse.cluster.default <- function(x, K1, chr = NULL, pos = NULL, ...) {
  # Produces a hierarchical clustering tree based on the input arrays.
  # Input: matrix K0, matrix K1
  # Output: matrix

  K0 <- x

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

  if(!all(is.integer(chr.idx),
          length(chr.idx) == nrow(K0),
          is.integer(pos),
          length(pos) == nrow(K0)))
    stop("Wrong input for fuse_cluster_R")

  # All set, calling the fuse_cluster_R now
  tree <- .Call('fuse_cluster_R', K0, K1, chr.idx, pos)
  tree <- .Call('sort_tree_R', tree)

  return(`colnames<-`(tree, c("m1", "m2", "logl_tot", "logl_merge", "genomic_dist")))
}


#' @export
fuse.cluster.matrix <- function(x, K1, chr = NULL, pos = NULL, ...) {

  if (missing(K1) || is.null(K1)) {
    stop("For matrix input, 'K1' must be supplied", call. = FALSE)
  }

  if (is.null(chr)) chr <- rep("chr1", nrow(x))
  if (is.null(pos)) pos <- seq_len(nrow(x))

  if (!all(dim(x) == dim(K1),
           length(chr) == nrow(x),
           length(pos) == nrow(x))) {
    stop("Incorrect input dimensions", call. = FALSE)
  }

  fuse.cluster.default(x, K1, chr, pos)
}


