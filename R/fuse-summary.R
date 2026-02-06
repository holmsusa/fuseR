#' Summarize FUSE Segmentation Results
#'
#' @description
#' Summarizes FUSE segmentation results into one row per segment,
#' including genomic coordinates, CpG count, segment length, average
#' methylation (beta), and stability flag based on likelihood testing.
#' Also returns per-sample methylation estimates for each segment.
#' Result can be visualized using plot(result).
#'
#' @param K0 Integer or numeric matrix of unmethylated counts.
#' @param K1 Integer or numeric matrix of methylated counts.
#' @param chr Character vector giving the chromosome for each site.
#' @param pos Numeric vector giving genomic coordinates for each site.
#' @param segments Integer vector giving segment IDs for each site in K0 and K1.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{summary}{A data frame with one row per segment and columns:}
#'     \describe{
#'       \item{Segment}{Segment ID}
#'       \item{Chr}{Chromosome}
#'       \item{Start}{Start genomic coordinate}
#'       \item{End}{End genomic coordinate}
#'       \item{CpGs}{Number of CpGs in the segment}
#'       \item{Length}{Genomic length (End - Start + 1)}
#'       \item{Beta}{Average methylation across samples and CpGs}
#'       \item{Coherent}{Logical indicator (TRUE if segment is coherently methylated, else FALSE)}
#'     }
#'   \item{betas_per_segment}{Matrix of per-sample methylation estimates for each segment
#'   (rows = segments, columns = samples).}
#'   \item{raw_beta}{Average beta per CpG site, used for plotting.}
#'   \item{raw_pos}{Genomic position for every given CpG site, used for plotting.}
#' }
#'
#' @examples
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
#' segments <- fuse.cut.tree(tree, 4)
#' res <- fuse.summary(K0, K1, rep("chr1", nrow(K0)), 1:nrow(K0), segments)
#' head(res$summary)
#' head(res$betas_per_segment)
#'
#'
#' @export
fuse.summary <- function(K0, K1, chr, pos, segments) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Summary support requires the 'stats' package.\n",
      "Install it with install.packages('stats').",
      call. = FALSE
    )
  }

  # --- Input validation ---
  if(!all((is.matrix(K0) || methods::is(K0, "DelayedArray")),
          (is.matrix(K1) || methods::is(K1, "DelayedArray")),
          all(dim(K0) == dim(K1))))
    stop("Wrong K0/K1 format")


  # --- Check colnames consistency ---
  if (!is.null(colnames(K0)) && !is.null(colnames(K1))) {
    if (!identical(colnames(K0), colnames(K1))) {
      stop("Column names of K0 and K1 must match exactly (including order).", call. = FALSE)
    }
  }

  # --- Initialization ---
  N <- ncol(K0)
  n_segments <- length(unique(segments))

  # --- Helper function ---
  safe_sum <- function(x) sum(x, na.rm = TRUE)

  # --- Model 1: Independent CpG model ---
  phat_ind <- .Call("bino_div_R", K1, K1 + K0)
  logL_ind <- -(.Call("bino_xlogy_R", K1, phat_ind) + .Call("bino_xlogy_R", K0, 1 - phat_ind))
  LL1 <- tapply(rowSums(logL_ind), segments, safe_sum)
  df1 <- tapply(rep(N, nrow(K0)), segments, safe_sum)

  # --- Model 2: Common segment model ---
  K0_seg <- rowsum(K0, group = segments, na.rm = TRUE)
  K1_seg <- rowsum(K1, group = segments, na.rm = TRUE)
  group_sizes <- as.numeric(table(segments))

  phat_seg <- .Call("bino_div_R", K1_seg, K1_seg + K0_seg)
  logL_seg <- -(.Call("bino_xlogy_R", K1_seg, phat_seg) + .Call("bino_xlogy_R", K0_seg, 1 - phat_seg))
  LL2 <- rowSums(logL_seg)
  df2 <- rep(N, n_segments)

  # --- Compute p-values and stability flag (using base R) ---
  pvals <- 1 - stats::pchisq(-2 * (LL1 - LL2), df = df1 - df2)
  coherent_flag <- pvals > 0.05

  # --- Compute average beta per segment ---
  total_K0 <- rowSums(K0_seg, na.rm = TRUE)
  total_K1 <- rowSums(K1_seg, na.rm = TRUE)
  beta_seg <- .Call("bino_div_R", total_K1, total_K0 + total_K1)

  # --- Compute genomic summaries ---
  segment_start <- tapply(pos, segments, min)
  segment_end   <- tapply(pos, segments, max)
  segment_chr   <- tapply(chr, segments, unique)
  segment_len   <- segment_end - segment_start + 1
  segment_cpgs  <- as.numeric(table(segments))

  # --- Combine into final summary data frame ---
  summary_df <- data.frame(
    Segment = paste(segment_chr, segment_start, sep = "."),
    Chr = as.character(unlist(segment_chr)),
    Start = as.integer(unlist(segment_start)),
    End = as.integer(unlist(segment_end)),
    CpGs = as.integer(segment_cpgs),
    Length = as.integer(segment_len),
    Beta = as.numeric(beta_seg),
    Coherent = coherent_flag
  )

  rownames(summary_df) <- NULL

  # --- Ensure betas_per_segment has correct dimnames ---
  if (!is.null(colnames(K1))) {
    colnames(phat_seg) <- colnames(K1)
  } else if (!is.null(colnames(K0))) {
    colnames(phat_seg) <- colnames(K0)
  }

  # Give rows the same identifiers as summary_df
  rownames(phat_seg) <- summary_df$Segment

  # --- Compute raw betas (per CpG) ---
  raw_beta <- rowSums(K1) / rowSums(K1 + K0)

  # --- Return both summary and per-segment betas ---
  out <- list(
    summary = summary_df,
    betas_per_segment = phat_seg,
    raw_beta = raw_beta,
    raw_pos = pos
  )

  class(out) <- "fuse_summary"
  return(out)

}





