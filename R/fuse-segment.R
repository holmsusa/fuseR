#' Full FUSE segmentation pipeline
#'
#' @description
#' Performs the full FUSE segmentation workflow:
#' hierarchical clustering, model selection, tree cutting,
#' and genomic segment summarization.
#'
#' @details
#' `fuse.segment()` is an S3 generic with methods for:
#' \describe{
#'   \item{\code{matrix}}{Raw count matrices (K0, K1) with genomic annotation.}
#'   \item{\code{BSseq}}{Bioconductor \code{BSseq} objects.}
#'   \item{\code{methrix}}{Bioconductor \code{methrix} objects (supports DelayedMatrix).}
#' }
#'
#' @param x Input object. One of:
#' \describe{
#'   \item{matrix}{Unmethylated count matrix (K0).}
#'   \item{BSseq}{A \code{BSseq} object.}
#'   \item{methrix}{A \code{methrix} object.}
#' }
#'
#' @param ... Additional arguments depending on input type:
#' \describe{
#'   \item{K1}{Methylated count matrix (required if \code{x} is a matrix).}
#'   \item{chr}{Chromosome labels, one per CpG (matrix input only).}
#'   \item{pos}{Genomic positions, one per CpG (matrix input only).}
#'   \item{method}{Information criterion for model selection:
#'     \code{"BIC"} (default) or \code{"AIC"}.}
#' }
#'
#' For internal use, `x` corresponds to the unmethylated count matrix (`K0`).
#'
#' @return
#' An object of class \code{fuse_summary}, containing:
#' \describe{
#'   \item{summary}{Data frame with one row per genomic segment.}
#'   \item{betas_per_segment}{Matrix of per-sample methylation estimates.}
#'   \item{raw_beta}{Per-CpG methylation estimates.}
#'   \item{raw_pos}{Genomic positions of CpGs.}
#' }
#'
#' @section Automatic data extraction:
#' For \code{BSseq} objects:
#' \itemize{
#'   \item Methylated counts are obtained via \code{getCoverage(x, "M")}
#'   \item Unmethylated counts via \code{getCoverage(x, "Cov") - M}
#'   \item Chromosome and position from \code{rowRanges(x)}
#' }
#'
#' For \code{methrix} objects:
#' \itemize{
#'   \item Methylated counts via \code{get_matrix(x, "M")}
#'   \item Total coverage via \code{get_matrix(x, "C")}
#'   \item Unmethylated counts computed as \code{C - M}
#'   \item Genomic coordinates extracted from locus metadata
#' }
#'
#'
#' @export
fuse.segment <- function(x, ...) {
  dots <- list(...)

  if (is.matrix(x) && length(dots) > 0 && is.null(names(dots)[1])) {
    dots$K1 <- dots[[1]]
    dots[[1]] <- NULL
  }

  UseMethod("fuse.segment", x)
}


#' @rdname fuse.segment
#' @param K1 Methylated count matrix (required if \code{x} is a matrix)
#' @param chr Chromosome labels, one per CpG (matrix input only)
#' @param pos Genomic positions, one per CpG (matrix input only)
#' @param method Information criterion for model selection: "BIC" (default) or "AIC"
#' @export
fuse.segment.default <- function(x, K1, chr, pos, method = c("BIC", "AIC"), ...) {
  blocks <- .materialize_by_chr(x, K1, chr)

  if(!all((is.numeric(pos) || is.integer(pos)),
          (is.character(chr)),
          length(pos) == length(chr)))
    stop("Wrong input format!")

  method <- match.arg(method)

  result <- list(summary = NULL, betas_per_segment = NULL, raw_beta = NULL, raw_pos = NULL)

  for (block in blocks) {
    K0_block <- block$K0
    K1_block <- block$K1

    tree <- fuse.cluster(K0_block, K1_block, chr[block$idx], pos[block$idx])
    tree[, 3] <- cumsum(tree[, 3])
    k_opt <- number.of.clusters(tree, ncol(K0_block), method)
    segments <- fuse.cut.tree(tree, k_opt)

    # Assert all elements have the same lengths
    stopifnot(length(chr[block$idx]) == nrow(K0_block))
    stopifnot(length(pos[block$idx]) == nrow(K0_block))
    stopifnot(length(segments) == nrow(K0_block))

    block_result <- fuse.summary(K0_block, K1_block, chr[block$idx], pos[block$idx], segments)

    result$summary <- rbind(result$summary, block_result$summary)
    result$betas_per_segment <- rbind(result$betas_per_segment, block_result$betas_per_segment)
    result$raw_beta <- c(result$raw_beta, block_result$raw_beta)
    result$raw_pos <- c(result$raw_pos, block_result$raw_pos)
  }

  class(result) <- "fuse_summary"

  result
}



#' @export
fuse.segment.matrix <- function(x, K1, chr = NULL, pos = NULL,
                                method = c("BIC", "AIC"), ...) {

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

  method <- match.arg(method)

  fuse.segment.default(x, K1, chr, pos, method)
}





