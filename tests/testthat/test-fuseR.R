test_that("fuse.cluster returns matrix with correct column names", {
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)

  result <- fuse.cluster(K0, K1)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 5)
  expect_equal(
    colnames(result),
    c("m1", "m2", "logl_tot", "logl_merge", "genomic_dist")
  )
})

test_that("fuse.sort.tree preserves matrix structure", {
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  tree <- fuse.cluster(K0, K1, sort = FALSE)

  sorted <- fuse.sort.tree(tree)

  expect_true(is.matrix(sorted))
  expect_equal(dim(sorted), dim(tree))
})

test_that("fuse.cutree returns integer vector of correct length", {
  tree <- matrix(c(
    -1, -2,  49.5,  49.5,  1.1,
    -3, -4,  78.5,  78.5,  1.1,
    -5, -6, 147.0, 147.0,  1.1,
    1,  2,  72.9, 201.0,  1.1,
    4,  3, 106.3, 454.4,  1.1
  ), ncol = 5, byrow = TRUE)

  segments <- fuse.cutree(tree, 3)

  expect_true(is.integer(segments))
  expect_equal(length(segments), 6)
})
