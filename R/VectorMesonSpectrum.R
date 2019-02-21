# Compute the spectrum of the rho meson
#' @export
VectorMesonSpectra <- function(nEigenvalues = 5)
{
  s <- computeSpectrum(z, urho, nEigen = nEigenvalues)
  s
}
