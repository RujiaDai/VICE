
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VICE

<!-- badges: start -->
<!-- badges: end -->

The goal of VICE is to calculate gene-level variability in single-cell or single-nuclei RNAseq data by constructing pseudo-replicates.

<img src="CV calculation.png" width="40%" />

## Installation

You can install the development version of VICE from
[GitHub](https://github.com/) with:

``` r
# # install.packages("devtools")
devtools::install_github("RujiaDai/VICE")
```

## Example

This is a basic example of using VICE:

``` r
library(VICE)
data(cmat)
data(cmeta)
cvlist <- get_cv_for_replicates(cmeta, cmat, 3)
cvplot(cvlist, "s1", "c1")
```

The parameter of VICE includes: the count matrix from sc/snRNAseq study `cmat` (gene by cell), the metadata of cells `cmeta` (cell by feature, "sample" and "celltype" must be provided), number of pseudo-replicates `k`.

The application of VICE on multiple sc/snRNAseq data can be found in <https://www.biorxiv.org/content/10.1101/2024.04.12.589216v1>.
