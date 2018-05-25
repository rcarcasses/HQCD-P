#' @export
ppDSigma <- function() DSigma('pp')

#' predict the values
#' @param dsigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param hs A data frame which contains the coefficients of  h(J) = h0 + h1 * (J-1) + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
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