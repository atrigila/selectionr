
<!-- README.md is generated from README.Rmd. Please edit that file -->

selectionr <img src="/Users/Usuario/Desktop/tests/hex.png" align="right" width="120" />
=======================================================================================

<!-- badges: start -->
<!-- badges: end -->

[![R build
status](https://github.com/tidyverse/ggplot2/workflows/R-CMD-check/badge.svg)](https://github.com/tidyverse/ggplot2/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tidyverse/ggplot2/master.svg)](https://codecov.io/github/tidyverse/ggplot2?branch=master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/ggplot2)](https://cran.r-project.org/package=ggplot2)

Overview
--------

`selectionr` is a personal package generated to ease repetitive taks in
positive selection screenings. The goal of selectionr is to simplify
general manipulation of multiple sequence alignments for positive
selection analysis within R.

Installation
------------

You can install the development version of selectionr from
[GitHub](https://github.com) with:

    devtools::install_github("atrigila/selectionr")

Example
-------

These are basic examples of the functionality of the package:

    library(selectionr)
    foxp2.sequences <- download.orthologues(gene.name = "FOXP2", target_species = "human", target_taxon = 9443)
