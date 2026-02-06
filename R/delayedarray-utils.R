.materialize_by_chr <- function(K0, K1, chr) {
  split_idx <- split(seq_len(nrow(K0)), chr)

  lapply(split_idx, function(idx) {
    list(
      K0 = as.matrix(K0[idx, , drop = FALSE]),
      K1 = as.matrix(K1[idx, , drop = FALSE]),
      idx = idx
    )
  })
}
