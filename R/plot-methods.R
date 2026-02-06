#' Plot method for FUSE segmentation results
#'
#' @description
#' Plotting method for fuse_summary.
#'
#' @details
#' Raw CpG-level methylation values are shown as grey points.
#' Segment-level methylation is shown as horizontal bars
#' (red = hypermethylated, blue = hypomethylated).
#'
#' @param x A fuse_summary object
#' @param ... Additional arguments
#' @param segments_to_plot Integer vector of segment indices
#'
#' @export
plot.fuse_summary <- function(x, ..., segments_to_plot = 1:50) {

  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop(
      "Plot support requires the 'graphics' package.\n",
      "Install it with install.packages('graphics').",
      call. = FALSE
    )
  }

  # If default
  if(identical(segments_to_plot, c(1:50))){
    segments_to_plot <- 1:min(50, nrow(x$summary))
  }

  if(!(all(segments_to_plot %in% 1:nrow(x$summary)))) {
    stop("segments_to_plot not in range")
  }


  summary <- x$summary[segments_to_plot, ]
  betas <- x$betas_per_segment[segments_to_plot, , drop = FALSE]

  betas_df <- data.frame(
    beta = rowMeans(betas),
    start = summary$Start,
    end = summary$End
  )

  # Indices of points to plot
  points_to_plot <- sum(x$summary$CpGs[1:(rev(segments_to_plot)[1]-1)]):sum(summary$CpGs)

  points_df <- data.frame(
    pos = x$raw_pos[points_to_plot],
    beta = x$raw_beta[points_to_plot]
  )

  plot(
    NA,
    xlim = range(c(betas_df$start, betas_df$end)),
    ylim = range(c(0,1)),
    xlab = "Genomic Position",
    ylab = "Beta",
    main = "Segments"
  )

  graphics::grid()
  graphics::abline(h = 0.5, col = "gray50", lty = 2, lwd = 1.5)

  # --- NEW: raw CpG beta points ---
  graphics::points(
    points_df$pos,
    points_df$beta,
    col = "grey",
    pch = 16

  )

  # Segment lines
  for (i in seq_len(nrow(summary))) {
    col <- if (betas_df$beta[i] > 0.5) "red" else "blue"
    graphics::segments(
      x0 = betas_df$start[i],
      x1 = betas_df$end[i],
      y0 = betas_df$beta[i],
      y1 = betas_df$beta[i],
      col = col,
      lwd = 4
    )
  }
}
