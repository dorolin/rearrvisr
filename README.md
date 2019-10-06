
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![Travis build status](https://travis-ci.org/dorolin/rearrvisr.svg?branch=master)](https://travis-ci.org/dorolin/rearrvisr) <!-- badges: end -->

rearrvisr
=========

The R package rearrvisr provides functions to identify and visualize inter- and intrachromosomal translocations and inversions between a focal genome and an ancestral genome reconstruction, or two extant genomes. Rearrangements, breakpoints, and synteny blocks are identified along the focal genome and output in a tabular format. Rearrangements and synteny blocks can be visualized along the focal genome by two graphical functions.

Installation
------------

You can install rearrvisr from GitHub using the devtools package:

``` r
## install and load devtools
install.packages("devtools")
library(devtools)

## install and load rearrvisr
install_github("dorolin/rearrvisr")
library(rearrvisr)
```

Example
-------

For a brief overview of the functionalities of the package, type

``` r
?rearrvisr
```

A more detailed description of the package functions and a tutorial is provided in the package vignette.
