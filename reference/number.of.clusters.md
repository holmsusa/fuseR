# Find Optimal Number of Clusters

Determines the optimal number of clusters to cut a hierarchical
clustering tree, based on the selected information criterion (e.g., BIC
or AIC).

## Usage

``` r
number.of.clusters(tree, n, method = c("BIC", "AIC"))
```

## Arguments

- tree:

  Clustering tree of class `hclust`.

- n:

  Number of samples in the original data.

- method:

  Information criterion method. One of `"BIC"` or `"AIC"`.

## Value

An integer representing the optimal number of clusters.

## Examples

``` r
# Example: Determine number of clusters in dummy data set
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
k <- number.of.clusters(tree, ncol(K0), 'BIC')
k
#> [1] 4
```
