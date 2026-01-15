# fishash

## Introduction

`fishash` is an R package for calling guides in perturbseq data from
UMI counts, based on treating the count matrix as a contingency table.

For each cell barcode and guide barcode, `fishash` tests how likely
the 2 barcodes are to co-occur in reads.  More specifically, it tests
whether the 2 barcodes have an odds ratio greater than 1 using a
one-sided Fisher's exact test. The method also includes a novel
procedure for correcting hidden confounding due to Simpson's paradox,
and performs a block-dependence-aware multiple testing correction
(assuming that tests from different cells are independent, but tests
within a cell are dependent).

## Installation

To install the package, by clone this repo, and run:
```{R}
devtools::install("/path/to/fishash")
```

## Basic usage

Given a count matrix `counts_mat` with guides for rows and cells for
columns, you can assign the guides calling the `fishash()` function:

```{R}
library(fishash)

# returns a SummarizedExperiment
res_fishash <- fishash(counts_mat)

# for the first few cells, print whether they received 0, 1, or 2+ guides:
head(colData(res_fishash)$demux_type)

# print the assigned guides for the first few cells
head(colData(res_fishash)$assignment)
```

For more options, see the help page:

```{R}
help(fishash)
```
