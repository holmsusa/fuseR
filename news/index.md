# Changelog

## fuseR 1.0.0 (2025-12-15)

### New features

- Added
  [`plot.fuse_summary()`](https://holmsusa.github.io/fuseR/reference/plot.fuse_summary.md)
  S3 method for visualization of FUSE segmentation results.
- Segment plots now display:
  - Horizontal segment-level methylation estimates
  - Per-segment average beta values
  - Optional per-CpG methylation values overlaid as background points
- Added support for flexible segment selection when plotting
  (`segments_to_plot`).

### Improvements

- Promoted package to first stable release (`1.0.0`).
- Improved input validation and error handling across exported
  functions.
- Improved S3 method registration and namespace imports.
- Improved vignette robustness and reproducibility during package
  checks.
- Updated documentation to reflect current plotting and summary
  capabilities.

### Bug fixes

- Fixed S3 dispatch issues affecting
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) during
  vignette building.
- Fixed lazy-load and installation issues caused by partial rebuilds.
- Fixed edge cases in beta computation when counts sum to zero.
- Fixed the swapped TRUE/FALSE labels for coherent segments in
  [`fuse.summary()`](https://holmsusa.github.io/fuseR/reference/fuse.summary.md).

------------------------------------------------------------------------

## fuseR 0.0.0.9000 (2025-11-20)

### Initial development release

- Implemented core FUSE pipeline:
  - [`fuse.cluster()`](https://holmsusa.github.io/fuseR/reference/fuse.cluster.md)
  - [`number.of.clusters()`](https://holmsusa.github.io/fuseR/reference/number.of.clusters.md)
  - [`fuse.cut.tree()`](https://holmsusa.github.io/fuseR/reference/fuse.cut.tree.md)
  - [`fuse.summary()`](https://holmsusa.github.io/fuseR/reference/fuse.summary.md)
  - [`fuse.segment()`](https://holmsusa.github.io/fuseR/reference/fuse.segment.md)
    (full pipeline wrapper)
- Included basic example data files for testing.
- Added introductory README and MIT license.
- Added full example workflow as an R Markdown vignette (`example.Rmd`).
- Enabled automatic pkgdown site deployment via GitHub Actions.
