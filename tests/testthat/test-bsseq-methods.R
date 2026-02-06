# --------- Works on bsseq objects  ----------------
test_that("fuse.segment works on BSseq objects", {
  skip_if_not_installed("bsseq")

  library(bsseq)

  set.seed(1)
  M <- matrix(sample(0:10, 200, TRUE), nrow = 50)
  Cov <- M + matrix(sample(0:10, 200, TRUE), nrow = 50)

  bs <- BSseq(
    chr = rep("chr1", 50),
    pos = seq_len(50),
    M = M,
    Cov = Cov,
    sampleNames = colnames(M)
  )

  res <- fuse.segment(bs)

  expect_s3_class(res, "fuse_summary")
})

