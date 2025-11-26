# ![fuseR hexagon](reference/figures/fuse_logo_meth.svg) fuseR

**fuseR** implements FUSE: **FUnctional SEgmentation of DNA methylation
data** through hierarchical clustering.

------------------------------------------------------------------------

## Features

- Hierarchical clustering based on methylation count matrices
- Detection of optimal number of clusters
- Summarization of segments, including stability flag
- Per-segment methylation estimates per sample

------------------------------------------------------------------------

## Installation

### From GitHub (Development Version)

``` r
# Install devtools or remotes if needed
install.packages("devtools")
# or
install.packages("remotes")

# Install fuseR from GitHub
devtools::install_github("holmsusa/fuseR")
# or
remotes::install_github("holmsusa/fuseR")
```

## System Requirements

- R version ≥ 4.0
- C++ toolchain for native code compilation

You may need platform-specific tools:

- macOS: Xcode Command Line Tools (xcode-select –install)
- Linux: build-essential
- Windows: Rtools

## Quick Example

``` r
library(fuseR)
set.seed(1234)

# Generate sample data
# Unmethylated counts, T's
K0 <- matrix(
  rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)
# Methylated counts, C's
K1 <- matrix(
  rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
  nrow = 100, byrow = TRUE
)

# Perform segmentation
segment_result <- fuse.segment(
  K0, K1, 
  chr = sub("\\..*$", "", rownames(K0)), 
  pos = as.numeric(sub("^.*\\.", "", rownames(K0)))
)

# Access summary and per-segment betas
head(segment_result$summary)
head(segment_result$betas_per_segment)
```

Check out a **full example workflow** in the
[vignette](https://holmsusa.github.io/fuseR/articles/example.html).

## License

This package is licensed under the MIT License. See
[LICENSE](https://holmsusa.github.io/fuseR/LICENSE) for details.

## Author

Susanna Holmstrom
