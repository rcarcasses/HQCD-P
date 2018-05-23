[![license](https://img.shields.io/github/license/rcarcasses/HQCD-P.svg)]()
[![GitHub language count](https://img.shields.io/github/languages/count/rcarcasses/HQCD-P.svg)]()
[![GitHub issues](https://img.shields.io/github/issues/rcarcasses/HQCD-P.svg)]()


# HQCD-P
An R package for the Pomeron in Holographic QCD.

This work is a collection of routines developed on the following papers:
- Soft Pomeron in Holographic QCD [arXiv:1508.00008](https://arxiv.org/abs/1508.00008)
- Unity of pomerons from gauge/string duality [arXiv:1704.08280](https://arxiv.org/abs/1704.08280)
- Non-minimal coupling contribution to DIS at low $x$ in Holographic QCD [arXiv:1804.07778](https://arxiv.org/abs/1804.07778)

# Dependencies

## Redis
This packages requires [redis](https://redis.io/) running on your computer or in a local accessible one. Redis is an *in-memory* database which is used to cache long computations. You can provide the IP of the redis instance in the `init()` function or disable it. However it is highly recommended that you use the cache system provided.

## `schrodinger` package
This package uses, among other packages which are available on CRAN, the [schrodinger](https://github.com/rcarcasses/schrodinger) package to find the eigenfunctions and eigenvalues of the associated schrodinger problem for the kernel of the pomeron. Make sure you have it installed it before using this. Make sure you have `armadillo` binaries in your system.

# Acknowledgements
This research received funding from the [European Union] 7th Framework Programme (Marie Curie Actions) under grant agreement 317089 (GATIS), from the grant CERN/FIS-NUC/0045/2015 and from the Simons Foundation grant 488637 (Simons collaboration on the Non-perturbative bootstrap). Centro de Física do Porto is partially funded by the Foundation for Science
and Technology of Portugal (FCT).
