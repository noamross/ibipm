[![](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/noamross/ibipm/blob/master/LICENSE.md)
[![](https://licensebuttons.net/l/by/3.0/80x15.png)](https://github.com/noamross/ibipm/blob/master/LICENSE.md)

[![DOI](https://zenodo.org/badge/6023/noamross/ibipm.svg)](https://zenodo.org/badge/latestdoi/6023/noamross/ibipm)

# Individual-Based Integral Projection Models

This is a repository contains the code for:

> Schreiber, Sebastian and Noam Ross (2015). Individual-Based Integral Projection
Models: The Role Of Size-Structure On Extinction Risk And Establishment
Success

`ibipm.Rnw` is an R-latex file of the manuscript.

`ibipm.pdf` is the compiled manuscript.

`ipm.R` contains the code for the IPM model described, and is sourced in
`ibipm.Rnw`.

`ibm-simulator.R` contains code for a full equivalent individual-based model simulator,
and is also sourced in `ibipm.Rnw`.

`ibipm.R` contains the code used in `ipipm.Rnw`, extracted using `purl("ibipm.Rnw")`.

`references.bib` contains the references used in the paper.
