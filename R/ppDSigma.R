#' @export
ppDSigma <- function() DSigma('pp')

#' @export
predict.ppDSigma <- function(dsigma, Izs, IzsBar, points, ...) {
  amplitude <- getAmplitude(dsigma, Izs, IzsBar, points, ...)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / (16 * pi^2)) * GEVMinus2ToNB * (1 / W^4) * abs(amplitude)^2
}

#' @export
IzN.ppDSigma <- function(pp, J, coeffs) {
  # Ansatz to compute the integral related to the proton
  sum(((J - 1)^(0:(length(coeff) - 1))) * coeff)
}

#' @export
IzNBar.ppDSigma <- IzN.ppDSigma
