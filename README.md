# diffeRenTES
<!---
[![Build Status](https://travis-ci.com/mbraccini/diffeRenTES.svg?branch=master)](https://travis-ci.com/mbraccini/diffeRenTES)
![Build Status GitHub Actions](https://github.com/mbraccini/diffeRenTES/actions/workflows/r.yml/badge.svg)
-->

<!-- badges: start -->
[![R-CMD-check](https://github.com/mbraccini/diffeRenTES/workflows/R-CMD-check/badge.svg)](https://github.com/mbraccini/diffeRenTES/actions) [![Codecov test coverage](https://codecov.io/gh/mbraccini/diffeRenTES/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mbraccini/diffeRenTES?branch=master)
<!-- badges: end -->


# What is it?
*diffeRenTES* (from French various, different, separate) is a package written in R for computing, starting from a Boolean network, the ATM (Attractor Transition Matrix) structure and tree-like structures based on the TES (Threshold Ergodic Sets)  concept that describe the process of differentiation.

TESs are the mathematical abstractions with which a powerful model of differentiation represents the different cell types arising during ontogenesis.
The aforementioned model of differentiation based on Boolean networks is firstly described in *"A Dynamical Model of Genetic Networks for Cell Differentiation" Villani M, Barbieri A, Serra R (2011) A Dynamical Model of Genetic Networks for Cell Differentiation. PLOS ONE 6(3): e17703. <https://doi.org/10.1371/journal.pone.0017703>*

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
