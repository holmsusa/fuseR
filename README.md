# fuseR

**fuseR** is an R package that implements FUSE: FUnctional SEgmentation of DNA methylation data through hierarchical clustering.

## Features
- Hierarchical clustering based on methylation count matrices
- Tree sorting for optimal ordering
- Tree segmentation based on a given number of clusters or information criteria

## Installation

### From GitHub (Development Version)

To install the latest version of `fuseR` from GitHub:

```r
# Install devtools if it's not already installed
install.packages("devtools")

# Then install fuseR from GitHub
devtools::install_github("holmsusa/fuseR")
```

## System Requirements
- R version ≥ 4.0
- C++ toolchain for native code compilation

You may need platform-specific tools: 
- macOS: Xcode Command Line Tools (xcode-select --install)
- Linux: build-essential, e.g., sudo apt install build-essential
- Windows: Rtools 

## Usage 
```r
library(fuseR)

# Generate sample data
K0 <- matrix(sample(1:200, 100, replace = TRUE), ncol = 5)
K1 <- matrix(sample(1:200, 100, replace = TRUE), ncol = 5)

# Perform clustering
tree <- fuse.cluster(K0, K1)

# Sort tree
sorted_tree <- fuse.sort.tree(tree)

# Cut tree into 3 clusters
segments <- fuse.cutree(sorted_tree, 3)
```

## License
This package is licensed under the MIT License. See the LICENSE file for details.

## Author
Susanna Holmstrom
