test_that("fuse.segment works on methrix objects", {

  skip_if_not_installed("methrix")
  skip_if_not_installed("GenomicRanges")

  data(methrix_data, package = "methrix")

  mx <- methrix_data

  res <- fuse.segment(mx)

  expect_type(res, "list")
  expect_true(all(c("summary", "betas_per_segment") %in% names(res)))

  expect_s3_class(res$summary, "data.frame")
  expect_true(is.matrix(res$betas_per_segment))

  expect_true(!is.null(attr(res, "k_opt")))
  expect_true(attr(res, "method") %in% c("BIC", "AIC"))
})

test_that("methrix genomic ranges are respected", {

  skip_if_not_installed("methrix")

  data(methrix_data, package = "methrix")

  loci <- methrix::get_matrix(methrix_data, type = "M", add_loci = TRUE)[, 1:3]
  res <- fuse.segment(methrix_data)

  expect_true(all(res$summary$start >= min(loci$start)))
  expect_true(all(res$summary$end   <= max(loci$start)))
})

