# Full FUSE Segmentation Pipeline

Performs the full FUSE segmentation workflow, from hierarchical
clustering to automatic determination of the optimal number of clusters
and summarization of segmentation results.

This function combines the steps:

1.  Hierarchical clustering via
    [`fuse.cluster`](https://holmsusa.github.io/fuseR/reference/fuse.cluster.md)

2.  Optimal number of segments determination via
    [`number.of.clusters`](https://holmsusa.github.io/fuseR/reference/number.of.clusters.md)

3.  Cutting the tree into clusters via
    [`fuse.cut.tree`](https://holmsusa.github.io/fuseR/reference/fuse.cut.tree.md)

4.  Summarizing the segmentation via
    [`fuse.summary`](https://holmsusa.github.io/fuseR/reference/fuse.summary.md)

## Usage

``` r
fuse.segment(K0, K1, chr, pos, method = c("BIC", "AIC"))
```

## Arguments

- K0:

  Integer or numeric matrix of unmethylated counts.

- K1:

  Integer or numeric matrix of methylated counts.

- chr:

  Character vector giving chromosome for each CpG site.

- pos:

  Numeric vector of genomic coordinates for each CpG site.

- method:

  Information criterion to use for determining optimal clusters. One of
  `"BIC"` or `"AIC"`. Defaults to `"BIC"`.

## Value

A list with two elements (the same structure as
[`fuse.summary`](https://holmsusa.github.io/fuseR/reference/fuse.summary.md)):

- summary:

  A data frame summarizing segments (chromosome, start/end, CpG count,
  etc.)

- betas_per_segment:

  Matrix of per-segment methylation estimates.

## Examples

``` r
# \donttest{
set.seed(1234)
K0 <- matrix(
  rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)
K1 <- matrix(
  rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)
chr <- rep("chr1", nrow(K0))
pos <- 1:nrow(K0)
res <- fuse.segment(K0, K1, chr, pos, method = "BIC")
head(res$summary)
#>   Segment  Chr Start End CpGs Length      Beta Coherent
#> 1  chr1.1 chr1     1  25   25     25 0.7523929     TRUE
#> 2 chr1.26 chr1    26  50   25     25 0.2355517     TRUE
#> 3 chr1.51 chr1    51  75   25     25 0.7523929     TRUE
#> 4 chr1.76 chr1    76 100   25     25 0.2355517     TRUE
# }
```
