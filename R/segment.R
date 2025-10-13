#' @useDynLib fuseR, .registration = TRUE
NULL

#' Full FUSE Segmentation Pipeline
#'
#' @description
#' Performs the full FUSE segmentation workflow, from hierarchical clustering
#' to automatic determination of the optimal number of clusters and summarization
#' of segmentation results.
#'
#' This function combines the steps:
#' \enumerate{
#'   \item Hierarchical clustering via \code{\link{fuse.cluster}}
#'   \item Optimal number of segments determination via \code{\link{number.of.clusters}}
#'   \item Cutting the tree into clusters via \code{\link{fuse.cut.tree}}
#'   \item Summarizing the segmentation via \code{\link{fuse.summary}}
#' }
#'
#' @param K0 Integer or numeric matrix of unmethylated counts.
#' @param K1 Integer or numeric matrix of methylated counts.
#' @param chr Character vector giving chromosome for each CpG site.
#' @param pos Numeric vector of genomic coordinates for each CpG site.
#' @param method Information criterion to use for determining optimal clusters.
#'   One of `"BIC"` or `"AIC"`. Defaults to `"BIC"`.
#' @param sort Logical; whether to sort the clustering tree optimally. Defaults to TRUE.
#'
#' @return
#' A list with two elements (the same structure as \code{\link{fuse.summary}}):
#' \describe{
#'   \item{summary}{A data frame summarizing segments (chromosome, start/end, CpG count, etc.)}
#'   \item{betas_per_segment}{Matrix of per-segment methylation estimates.}
#' }
#'
#' @examples
#' \donttest{
#' K0 <- matrix(sample(1:200, 125, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(1:200, 125, replace = TRUE), ncol = 5)
#' chr <- rep("chr1", nrow(K0))
#' pos <- 1:nrow(K0)
#'
#' res <- fuse_segment(K0, K1, chr, pos, method = "BIC")
#' head(res$summary)
#' }
#'
#' @export
fuse.segment <- function(K0, K1, chr, pos, method = c("BIC", "AIC"), sort = TRUE) {
  # --- Input validation ---
  stopifnot(
    is.matrix(K0),
    is.matrix(K1),
    all(dim(K0) == dim(K1)),
    is.character(chr),
    length(chr) == nrow(K0),
    is.numeric(pos),
    length(pos) == nrow(K0)
  )

  method <- match.arg(method)

  # --- Step 1: Hierarchical clustering ---
  tree <- fuse.cluster(K0, K1, sort = sort)

  # --- Step 2: Determine optimal number of clusters ---
  # Need the total likelihood on each step, so the sum of the changes in total likelihood
  tree[,3] <- cumsum(tree[,3])
  k_opt <- number.of.clusters(tree, ncol(K0), method = method)

  # --- Step 3: Cut tree into segments ---
  segments <- fuse.cut.tree(tree, k_opt)

  # --- Step 4: Summarize results ---
  result <- fuse.summary(K0, K1, chr, pos, segments)

  # Attach metadata for convenience
  attr(result, "k_opt") <- k_opt
  attr(result, "method") <- method
  attr(result, "sorted") <- sort

  return(result)
}
