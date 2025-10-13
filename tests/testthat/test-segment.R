# fuse.segment
test_that("fuse.segment runs end-to-end and returns valid structure", {
  skip_if_not(exists("fuse.cluster"))
  skip_if_not(exists("number.of.clusters"))
  skip_if_not(exists("fuse.cut.tree"))
  skip_if_not(exists("fuse.summary"))

  set.seed(123)
  K0 <- matrix(sample(1:20, 30, replace = TRUE), ncol = 3)
  K1 <- matrix(sample(1:20, 30, replace = TRUE), ncol = 3)
  chr <- rep("chr1", nrow(K0))
  pos <- 1:nrow(K0)

  res <- fuse.segment(K0, K1, chr, pos, method = "AIC")

  expect_type(res, "list")
  expect_named(res, c("summary", "betas_per_segment"))
  expect_true(is.data.frame(res$summary))
  expect_true(is.matrix(res$betas_per_segment))
  expect_true(!is.null(attr(res, "k_opt")))
  expect_true(attr(res, "k_opt") > 0)
})

test_that("fuse.segment handles invalid inputs gracefully", {
  K0 <- matrix(1:10, ncol = 2)
  K1 <- matrix(1:10, ncol = 2)
  chr <- rep("chr1", nrow(K0))
  pos <- 1:nrow(K0)

  # Dimension mismatch
  expect_error(fuse.segment(K0, matrix(1:12, ncol = 2), chr, pos))

  # Invalid chr/pos lengths
  expect_error(fuse.segment(K0, K1, chr[-1], pos))
  expect_error(fuse.segment(K0, K1, chr, pos[-1]))

  # Invalid method
  expect_error(fuse.segment(K0, K1, chr, pos, method = "XYZ"))
})






