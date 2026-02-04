#' Full FUSE segmentation pipeline
#'
#' @description
#' Performs the full FUSE segmentation workflow, including clustering,
#' model selection, tree cutting, and segment summarization.
#'
#' @param x Input data.
#'   If a matrix, interpreted as the unmethylated count matrix (K0).
#'   If a BSseq object, methylation counts are extracted automatically.
#'
#' @param ... Additional arguments depending on the input type:
#'   \describe{
#'     \item{K1}{Methylated count matrix (required if \code{x} is a matrix).}
#'     \item{chr}{Character vector of chromosome labels, one per CpG site.}
#'     \item{pos}{Numeric vector of genomic positions, one per CpG site.}
#'     \item{method}{Information criterion for model selection:
#'       either \code{"BIC"} (default) or \code{"AIC"}.}
#'   }
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{summary}{Data frame summarizing genomic segments.}
#'   \item{betas_per_segment}{Matrix of per-segment methylation estimates.}
#' }
#'
#' @export
fuse.segment <- function(x, ...) {
  UseMethod("fuse.segment")
}


#' @export
fuse.segment.default <- function(x, K1, chr, pos, method = c("BIC", "AIC"), ...) {
  blocks <- .materialize_for_fuse(x, K1)

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
    result$raw_beta <- rbind(result$raw_beta, block_result$raw_beta)
    result$raw_pos <- rbind(result$raw_pos, block_result$raw_pos)
  }

  attr(result, "k_opt") <- k_opt
  attr(result, "method") <- method
  class(result) <- "fuse_summary"

  result
}



#' @export
fuse.segment.matrix <- function(x, K1, chr, pos, ...) {
  if (missing(chr) || missing(pos)) {
    stop("For matrix input, use fuse.segment(K0, K1, chr, pos)", call. = FALSE)
  }
  fuse.segment.default(x, K1, chr, pos, ...)
}



