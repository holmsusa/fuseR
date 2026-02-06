# Full FUSE segmentation pipeline

Performs the full FUSE segmentation workflow: hierarchical clustering,
model selection, tree cutting, and genomic segment summarization.

## Usage

``` r
fuse.segment(x, ...)

# Default S3 method
fuse.segment(x, K1, chr, pos, method = c("BIC", "AIC"), ...)

# S3 method for class 'BSseq'
fuse.segment(x, method = c("BIC", "AIC"), ...)

# S3 method for class 'methrix'
fuse.segment(x, method = c("BIC", "AIC"), ...)
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

  Additional arguments depending on input type:

  K1

  :   Methylated count matrix (required if `x` is a matrix).

  chr

  :   Chromosome labels, one per CpG (matrix input only).

  pos

  :   Genomic positions, one per CpG (matrix input only).

  method

  :   Information criterion for model selection: `"BIC"` (default) or
      `"AIC"`.

  For internal use, `x` corresponds to the unmethylated count matrix
  (`K0`).

- K1:

  Methylated count matrix (required if `x` is a matrix)

- chr:

  Chromosome labels, one per CpG (matrix input only)

- pos:

  Genomic positions, one per CpG (matrix input only)

- method:

  Information criterion for model selection: "BIC" (default) or "AIC"

## Value

An object of class `fuse_summary`, containing:

- summary:

  Data frame with one row per genomic segment.

- betas_per_segment:

  Matrix of per-sample methylation estimates.

- raw_beta:

  Per-CpG methylation estimates.

- raw_pos:

  Genomic positions of CpGs.

## Details

`fuse.segment()` is an S3 generic with methods for:

- `matrix`:

  Raw count matrices (K0, K1) with genomic annotation.

- `BSseq`:

  Bioconductor `BSseq` objects.

- `methrix`:

  Bioconductor `methrix` objects (supports DelayedMatrix).

## Automatic data extraction

For `BSseq` objects:

- Methylated counts are obtained via `getCoverage(x, "M")`

- Unmethylated counts via `getCoverage(x, "Cov") - M`

- Chromosome and position from `rowRanges(x)`

For `methrix` objects:

- Methylated counts via `get_matrix(x, "M")`

- Total coverage via `get_matrix(x, "C")`

- Unmethylated counts computed as `C - M`

- Genomic coordinates extracted from locus metadata
