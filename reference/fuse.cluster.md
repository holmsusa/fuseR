# Perform Hierarchical Clustering on Methylation Data

Produces a hierarchical clustering tree based on the input matrices of
counts.

## Usage

``` r
fuse.cluster(x, ...)

# Default S3 method
fuse.cluster(x, K1, chr = NULL, pos = NULL, ...)

# S3 method for class 'BSseq'
fuse.cluster(x, ...)

# S3 method for class 'methrix'
fuse.cluster(x, ...)
```

## Arguments

- x:

  Input object. One of:

  matrix

  :   Unmethylated count matrix (K0).

  BSseq

  :   A `BSseq` object.

  methrix

  :   A `methrix` object.

- ...:

  Additional arguments if K0 is a matrix.

- K1:

  Methylated count matrix

- chr:

  Optional chromosome vector

- pos:

  Optional position vector

## Value

A clustering tree of class `hclust`.

## Examples

``` r
# Example: Clustering generated data
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
tree
#> 
#> Call:
#> "fuse.cluster(k0, k1, pos, chr)"
#> 
#> Cluster method   : fuse 
#> Distance         : fuse 
#> Number of objects: 100 
#> 
```
