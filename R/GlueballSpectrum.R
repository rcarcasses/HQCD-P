# Compute the spectrum of the scalar glueballs 0^{++}
#' @export
ScalarGlueballSpectra <- function(nEigenvalues = 5)
{
  s <- computeSpectrum(z, u0, nEigen = nEigenvalues)
  s
}
# Compute the spectrum of the scalar glueballs 2^{++}
#' @export
TensorGlueballSpectra <- function(nEigenvalues = 5)
{
  s <- computeSpectrum(z, u2, nEigen = nEigenvalues)
  s
}
