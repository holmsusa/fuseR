.get_fuse_inputs <- function(x, ...) {
  UseMethod(".get_fuse_inputs")
}

.get_fuse_inputs.BSseq <- function(x, ...) {
  stopifnot(requireNamespace("bsseq", quietly = TRUE))

  K1 <- bsseq::getCoverage(x, type = "M")
  K0 <- bsseq::getCoverage(x, type = "Cov") - K1

  gr <- SummarizedExperiment::rowRanges(x)

  list(
    K0 = K0,
    K1 = K1,
    chr = as.character(GenomicRanges::seqnames(gr)),
    pos = GenomicRanges::start(gr)
  )
}

.get_fuse_inputs.methrix <- function(x, ...) {
  stopifnot(requireNamespace("methrix", quietly = TRUE))
  stopifnot(requireNamespace("DelayedArray", quietly = TRUE))

  # get_matrix() returns DelayedMatrix if add_loci=FALSE
  K1 <- methrix::get_matrix(x, type = "M", add_loci = FALSE)
  K  <- methrix::get_matrix(x, type = "C", add_loci = FALSE)

  # Compute unmethylated counts as DelayedMatrix
  K0 <- K - K1

  loci <- methrix::get_matrix(x, type = "M", add_loci = TRUE)[, 1:3]


  # Keep only rows present in K0
  valid_idx <- seq_len(nrow(K0))

  list(
    K0  = K0,
    K1  = K1,
    chr = as.character(loci$chr),
    pos = as.integer(loci$start)
  )
}


#' @export
fuse.segment.BSseq <- function(x, method = c("BIC", "AIC"), ...) {
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop(
      "BSseq support requires the 'bsseq' package.\n",
      "Install it with BiocManager::install('bsseq')",
      call. = FALSE
    )
  }

  inputs <- .get_fuse_inputs(x)

  fuse.segment.default(
    x   = inputs$K0,
    K1  = inputs$K1,
    chr = inputs$chr,
    pos = inputs$pos,
    method = method
  )
}

#' @export
fuse.segment.methrix <- function(x, method = c("BIC", "AIC"), ...) {
  if (!requireNamespace("methrix", quietly = TRUE)) {
    stop(
      "methrix support requires the 'methrix' package.\n",
      "Install it with BiocManager::install('methrix')",
      call. = FALSE
    )
  }

  inputs <- .get_fuse_inputs(x)

  fuse.segment.default(
    x   = inputs$K0,
    K1  = inputs$K1,
    chr = inputs$chr,
    pos = inputs$pos,
    method = method
  )
}


