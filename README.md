# diffeRenTES
<!---
[![Build Status](https://travis-ci.com/mbraccini/diffeRenTES.svg?branch=master)](https://travis-ci.com/mbraccini/diffeRenTES)
![Build Status GitHub Actions](https://github.com/mbraccini/diffeRenTES/actions/workflows/r.yml/badge.svg)
-->

<!-- badges: start -->
[![R-CMD-check](https://github.com/mbraccini/diffeRenTES/workflows/R-CMD-check/badge.svg)](https://github.com/mbraccini/diffeRenTES/actions) [![Codecov test coverage](https://codecov.io/gh/mbraccini/diffeRenTES/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mbraccini/diffeRenTES?branch=master)[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)[![Made with](https://img.shields.io/badge/Made%20with-♥%20and%20R-%3CCOLOR%3E)](https://github.com/mbraccini/diffeRenTES/)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/diffeRenTES)](https://cran.r-project.org/package=diffeRenTES)[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/diffeRenTES)](https://cranlogs.r-pkg.org/badges/grand-total/diffeRenTES)
<!-- badges: end -->


# What is it?
*diffeRenTES* (from French various, different) is a package written in R that, starting from the Boolean networks of the *BoolNet* package, computes the ATM (Attractor Transition Matrix) structure and the tree-like structure, based on the TES concept, describing the cell differentiation process.

TESs (Threshold Ergodic Sets) are the mathematical abstractions that represent the different cell types arising during ontogenesis.

TESs and the powerful model of biological differentiation based on Boolean networks to which it belongs have been firstly described in *"A Dynamical Model of Genetic Networks for Cell Differentiation" Villani M, Barbieri A, Serra R (2011) A Dynamical Model of Genetic Networks for Cell Differentiation. PLOS ONE 6 (3): e17703. <https://doi.org/10.1371/journal.pone.0017703>*

# Package installation

This package requires the following R libraries:

- BoolNet
- igraph
- DOT

```r
install.packages("diffeRenTES_RELEASE_VERSION.tar.gz", repos = NULL, type="source")
```
Or, alternatively, directly from GitHub:

```r
library(devtools)
devtools::install_github("mbraccini/diffeRenTES")
```

# Quick tutorial

```r
#Boolean network generation by means of 'BoolNet' package
net <- BoolNet::generateRandomNKNetwork(10, 2)

#Attractors computation
attractors <- BoolNet::getAttractors(net) 

#Attractors Transition Matrix computation
ATM <- getATM(net, attractors)

#TESs computation
TESs <- getTESs(ATM)

#Saving the image of the TES-based differentiation tree
saveDifferentiationTreeToFile(TESs, "example.svg") 
```
