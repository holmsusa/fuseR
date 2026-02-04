# Summarize FUSE Segmentation Results

Summarizes FUSE segmentation results into one row per segment, including
genomic coordinates, CpG count, segment length, average methylation
(beta), and stability flag based on likelihood testing. Also returns
per-sample methylation estimates for each segment. Result can be
visualized using plot(result).

## Usage

``` r
fuse.summary(K0, K1, chr, pos, segments)
```

## Arguments

- K0:

  Integer or numeric matrix of unmethylated counts.

- K1:

  Integer or numeric matrix of methylated counts.

- chr:

  Character vector giving the chromosome for each site.

- pos:

  Numeric vector giving genomic coordinates for each site.

- segments:

  Integer vector giving segment IDs for each site in K0 and K1.

## Value

A list with two elements:

- summary:

  A data frame with one row per segment and columns:

- betas_per_segment:

  Matrix of per-sample methylation estimates for each segment (rows =
  segments, columns = samples).

- raw_beta:

  Average beta per CpG site, used for plotting.

- raw_pos:

  Genomic position for every given CpG site, used for plotting.

## Examples

``` r
set.seed(1234)
K0 <- matrix(
  rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)
K1 <- matrix(
  rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)
tree <- fuse.cluster(K0, K1)
segments <- fuse.cut.tree(tree, 4)
res <- fuse.summary(K0, K1, rep("chr1", nrow(K0)), 1:nrow(K0), segments)
head(res$summary)
#>   Segment  Chr Start End CpGs Length      Beta Coherent
#> 1  chr1.1 chr1     1  25   25     25 0.7523929    FALSE
#> 2 chr1.26 chr1    26  50   25     25 0.2355517    FALSE
#> 3 chr1.51 chr1    51  75   25     25 0.7523929    FALSE
#> 4 chr1.76 chr1    76 100   25     25 0.2355517    FALSE
head(res$betas_per_segment)
#>              [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> chr1.1  0.7895277 0.7677083 0.7163462 0.7583893 0.7068311 0.7875920 0.7548919
#> chr1.26 0.2435233 0.2213667 0.2137791 0.2308478 0.2693467 0.2171665 0.2216700
#> chr1.51 0.7895277 0.7677083 0.7163462 0.7583893 0.7068311 0.7875920 0.7548919
#> chr1.76 0.2435233 0.2213667 0.2137791 0.2308478 0.2693467 0.2171665 0.2216700
#>              [,8]
#> chr1.1  0.7444562
#> chr1.26 0.2651515
#> chr1.51 0.7444562
#> chr1.76 0.2651515

```
