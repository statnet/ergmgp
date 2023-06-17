# `ergmgp`: Tools for Modeling ERGM Generating Processes

`ergmgp` contains tools for modeling continuous time graph processes with equilibrium distributions in exponential family random graph model (ERGM) form (*ERGM generating processes*, or *EGPs*).  The `ergmgp` tools allow models to be specified in terms of their equilibria, using [`ergm`](https://github.com/statnet/ergm) formulas, together with additional information on the dynamic process.  A number of different processes are supported, including longitudinal ERGMs, continuum limits of temporal and separable temporal ERGMs, and competing rate stochastic actor-oriented model processes, as well as differential stability and change inhibition processes.

An overview of supported EGPs and pointers to other help pages can be obtained after loading the package with `help(ergmgp)`.  Additional information on EGPs can also be found at the reference below.

## Installing from CRAN

The easiest way to install the package is to use CRAN.  From within `R`, simply use

```
install.packages("ergmgp")
```

which will install `ergmgp` and its dependencies.  Calling `library(ergmgp)` will subsequently load the package, and away you go.

## Installing Directly from GitHub

To install from GitHub, first ensure that you have the `devtools` package installed and loaded. Then, type the following: 

```
install_github("statnet/ergmgp")
```
Alternately, cloning this repository and building/installing the package locally is another option. 

## References

Butts, Carter T.  (2023).  ``Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.'' *Journal of Mathematical Sociology*.  DOI: 10.1080/0022250X.2023.2180001  [[arXiv version](https://arxiv.org/abs/2203.06948)]
