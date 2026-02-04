# fuse.summary
# ------------------ Returns correct structure and dimensions ----------------------
test_that("fuse.summary returns correct structure and dimensions", {
  # --- Synthetic data setup ---
  set.seed(123)
  K0 <- matrix(sample(0:50, 20, replace = TRUE), ncol = 4)
  K1 <- matrix(sample(0:50, 20, replace = TRUE), ncol = 4)
  chr <- rep("chr1", nrow(K0))
  pos <- seq_len(nrow(K0))
  segments <- rep(1:4, each = 5)[1:nrow(K0)]  # 4 segments

  # --- Run function ---
  res <- fuse.summary(K0, K1, chr, pos, segments)

  # --- Basic structure ---
  expect_type(res, "list")
  expect_named(res, c("summary", "betas_per_segment", "raw_beta", "raw_pos"))

  # --- Summary structure ---
  sm <- res$summary
  expect_s3_class(sm, "data.frame")
  expect_true(all(c("Segment", "Chr", "Start", "End",
                    "CpGs", "Length", "Beta", "Coherent") %in% names(sm)))

  # One row per segment
  expect_equal(nrow(sm), length(unique(segments)))

  # Betas per segment: rows = segments, columns = samples
  expect_equal(dim(res$betas_per_segment),
               c(length(unique(segments)), ncol(K0)))

  # Types and ranges
  expect_true(all(is.numeric(sm$Beta)))
  expect_true(all(is.logical(sm$Coherent)))
  expect_true(all(sm$CpGs > 0))
  expect_true(all(sm$Length >= 1))
})

# ------------------ Consistent segment ranges  ----------------------
test_that("fuse.summary computes consistent segment ranges", {
  # --- Simple deterministic data ---
  K0 <- matrix(10, nrow = 6, ncol = 3)
  K1 <- matrix(5, nrow = 6, ncol = 3)
  chr <- rep("chr2", 6)
  pos <- c(1, 2, 3, 10, 11, 12)
  segments <- c(1, 1, 1, 2, 2, 2)

  res <- fuse.summary(K0, K1, chr, pos, segments)
  sm <- res$summary

  # Start/End must match min/max of pos per segment
  expect_equal(sm$Start, c(1, 10))
  expect_equal(sm$End,   c(3, 12))
  expect_equal(sm$Length, sm$End - sm$Start + 1)
  expect_equal(sm$CpGs, c(3, 3))
  expect_true(all(sm$Chr == "chr2"))

  # Beta roughly 5/(10+5) = 0.333...
  expect_true(all(abs(sm$Beta - 1/3) < 1e-8))
})

# ------------------ Handles invalid inputs ----------------------
test_that("fuse.summary handles invalid inputs", {
  K0 <- matrix(1:6, ncol = 2)
  K1 <- matrix(1:6, ncol = 2)
  chr <- rep("chr1", nrow(K0))
  pos <- 1:nrow(K0)
  segments <- rep(1, nrow(K0))

  # Non-matrix inputs
  expect_error(
    fuse.summary(as.data.frame(K0), K1, chr, pos, segments)
  )
  expect_error(
    fuse.summary(K0, as.data.frame(K1), chr, pos, segments)
  )

  # Dimension mismatch
  expect_error(fuse.summary(K0, matrix(1:8, ncol = 2), chr, pos, segments))

  # Invalid chr length
  expect_error(fuse.summary(K0, K1, chr[-1], pos, segments))

  # Invalid pos type or length
  expect_error(fuse.summary(K0, K1, chr, as.character(pos), segments))
  expect_error(fuse.summary(K0, K1, chr, pos[-1], segments))

  # Invalid segments type or length
  expect_error(fuse.summary(K0, K1, chr, pos, segments[-1]))
})

# ------------------ Stability flag local and reproducible ----------------------
test_that("fuse.summary stability flag is logical and reproducible", {
  set.seed(42)
  K0 <- matrix(sample(0:5, 16, replace = TRUE), nrow = 16, ncol = 4)
  K1 <- matrix(sample(0:5, 16, replace = TRUE), nrow = 16, ncol = 4)
  chr <- rep("chr3", nrow(K0))
  pos <- seq_len(nrow(K0))
  segments <- rep(1:4, each = 4)

  res1 <- fuse.summary(K0, K1, chr, pos, segments)
  res2 <- fuse.summary(K0, K1, chr, pos, segments)

  expect_equal(res1$summary$Coherent, res2$summary$Coherent)
  expect_type(res1$summary$Coherent, "logical")
})

