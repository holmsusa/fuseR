# Cut Hierarchical Clustering Tree into Clusters

Divides the clustering tree into a specified number of clusters.

## Usage

``` r
fuse.cut.tree(tree, k)
```

## Arguments

- tree:

  Clustering tree of class `hclust`.

- k:

  Number of clusters

## Value

A vector indicating which cluster each element in the original data
frame belonged to

## Examples

``` r
# Example: Cutting small tree in 2 segments
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

segments <- fuse.cut.tree(tree, 2)
segments
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [75] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
```
