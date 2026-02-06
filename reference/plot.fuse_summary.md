# Plot method for FUSE segmentation results

Plotting method for fuse_summary.

## Usage

``` r
# S3 method for class 'fuse_summary'
plot(x, ..., segments_to_plot = 1:50)
```

## Arguments

- x:

  A fuse_summary object

- ...:

  Additional arguments

- segments_to_plot:

  Integer vector of segment indices

## Details

Raw CpG-level methylation values are shown as grey points. Segment-level
methylation is shown as horizontal bars (red = hypermethylated, blue =
hypomethylated).
