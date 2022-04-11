# diffeRenTES
<!---
[![Build Status](https://travis-ci.com/mbraccini/diffeRenTES.svg?branch=master)](https://travis-ci.com/mbraccini/diffeRenTES)
![Build Status GitHub Actions](https://github.com/mbraccini/diffeRenTES/actions/workflows/r.yml/badge.svg)
-->

<!-- badges: start -->
[![R-CMD-check](https://github.com/mbraccini/diffeRenTES/workflows/R-CMD-check/badge.svg)](https://github.com/mbraccini/diffeRenTES/actions) [![Codecov test coverage](https://codecov.io/gh/mbraccini/diffeRenTES/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mbraccini/diffeRenTES?branch=master)[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)[![Made with](https://img.shields.io/badge/Made%20with-♥%20and%20R-%3CCOLOR%3E)](https://github.com/mbraccini/diffeRenTES/)
<!-- badges: end -->


# What is it?
*diffeRenTES* (from French various, different, separate) is a package written in R that, starting from the Boolean networks of the *BoolNet* package, computes the ATM (Attractor Transition Matrix) structure and the tree-like structure, based on the TES concept, describing the cell differentiation process.

TESs (Threshold Ergodic Sets) are the mathematical abstractions that represent the different cell types arising during ontogenesis.

TESs and the powerful model of biological differentiation based on Boolean networks to which it belongs have been firstly described in *"A Dynamical Model of Genetic Networks for Cell Differentiation" Villani M, Barbieri A, Serra R (2011) A Dynamical Model of Genetic Networks for Cell Differentiation. PLOS ONE 6 (3): e17703. <https://doi.org/10.1371/journal.pone.0017703>*

# Package installation

This package requires the folowing R libraries:

- BoolNet
- igraph
- DOT

```r
install.packages("diffeRenTES_RELEASE_VERSION.tar.gz", repos = NULL, type="source")
```


# Quick tutorial

```r
  	net <- BoolNet::generateRandomNKNetwork(10, 2)
	attractors <- BoolNet::getAttractors(net) 
	ATM <- computeATM(net, attractors) # attractors transition matrix computation
	TESs <- computeTESs(ATM) #TESs computation
	saveDotFileDifferentiationTree(TESs, "exampleTree") #saving the TES-based differentiation tree
```
