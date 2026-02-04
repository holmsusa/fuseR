# Full FUSE segmentation pipeline

Performs the full FUSE segmentation workflow, including clustering,
model selection, tree cutting, and segment summarization.

## Usage

``` r
fuse.segment(x, ...)
```

## Arguments

- x:

  Input data. If a matrix, interpreted as the unmethylated count matrix
  (K0). If a BSseq object, methylation counts are extracted
  automatically.

- ...:

  Additional arguments depending on the input type:

  K1

  :   Methylated count matrix (required if `x` is a matrix).

  chr

  :   Character vector of chromosome labels, one per CpG site.

  pos

  :   Numeric vector of genomic positions, one per CpG site.

  method

  :   Information criterion for model selection: either `"BIC"`
      (default) or `"AIC"`.

## Value

A list with two elements:

- summary:

  Data frame summarizing genomic segments.

- betas_per_segment:

  Matrix of per-segment methylation estimates.
