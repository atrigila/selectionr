
<!-- README.md is generated from README.Rmd. Please edit that file -->

selectionr
==========

<!-- badges: start -->
<!-- badges: end -->

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
