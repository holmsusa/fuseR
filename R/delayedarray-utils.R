.materialize_for_fuse <- function(K0, K1, max_rows = 1e6) {
  if (!methods::is(K0, "DelayedArray")) {
    return(list(list(K0 = as.matrix(K0), K1 = as.matrix(K1), idx = seq_len(nrow(K0)))))
  }

  # Block-wise realization
  blocks <- beachmat::rowBlockApply(
    K0,
    FUN = function(idx) {
      list(
        K0 = as.matrix(K0[idx, , drop = FALSE]),
        K1 = as.matrix(K1[idx, , drop = FALSE]),
        idx = idx
      )
    }
  )

  blocks
}
