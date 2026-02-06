# <img src="man/figures/fuse_logo.svg" alt="fuseR hexagon" align="right" height="180" style="margin-left: 0.5em" /> fuseR

**fuseR** implements FUSE: **FUnctional SEgmentation of DNA methylation data** through hierarchical clustering.

---


## Features

- Hierarchical clustering based on methylation count matrices
- Detection of optimal number of clusters
- Summarization of segments, including stability flag
- Per-segment methylation estimates per sample

---

## Installation

### From GitHub (Development Version)
Either using remotes: (recommended)

```r
# Install remotes if needed
install.packages("remotes")

# Install fuseR from GitHub
remotes::install_github("holmsusa/fuseR")
```

or using devtools:

```r
# Install devtools if needed
install.packages("devtools")

# Install fuseR from GitHub
devtools::install_github("holmsusa/fuseR")
```

## System Requirements
- R version ≥ 4.0
- C++ toolchain for native code compilation

You may need platform-specific tools: 

- macOS: Xcode Command Line Tools (xcode-select --install)
- Linux: build-essential
- Windows: Rtools 

## Supported input types

`fuse.segment()` supports the following input formats:

- **Matrix input**
  - Unmethylated counts (`K0`)
  - Methylated counts (`K1`)
  - Explicit chromosome (`chr`) and position (`pos`) vectors

- **BSseq objects** (Bioconductor)
  - Counts and genomic coordinates are extracted automatically
  - Requires **bsseq** package

- **methrix objects** (Bioconductor)
  - Supports large datasets via `DelayedMatrix`
  - Counts and genomic coordinates are extracted automatically
  - Requires **methrix** and **DelayedArray** packages
  
Install needed packages with 
```r
BiocManager::install(c("bsseq", "methrix". "DelayedArray"))
```


## Quick Example 
```r
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

Check out a **full example workflow** in the [vignette](https://holmsusa.github.io/fuseR/articles/example.html).


## License
This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Author
Susanna Holmstrom
