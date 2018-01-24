# HQCD-P
A package for the Pomeron in Holographic QCD

# Depencences

## Redis dependence
This packages requires [redis](https://redis.io/) running on your computer or in a local accessible one. Redis is an *in memory* database which is used to cache long computations. You can provide the IP of the redis instance in the `init()` function or disable it. However it is highly recommended that you use the cache system provided.

## schrodinger package
This package uses, among other packages which are available on CRAN, the [schrodinger](https://github.com/rcarcasses/schrodinger) package to find the eigenfunctions and eigenvalues of the associated schrodinger problem for the kernel of the pomeron. Make sure you have it installed it before using this. Make sure you have `armadillo` binaries in your system.
