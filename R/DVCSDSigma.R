#' @export
DVCSDSigma <- function() DSigma('DVCSDSigma')

getExternalStateFactor.DVCSDSigma <- function(dvcs, Q2 = Q2, ...) {
  getU1NNMode(Q2 = Q2)$fQ
}

#' predict the values
#' @param dsigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.DVCSDSigma <- function(dsigma, fns, gs, points, ...) {
  # compute the amplitude
  amplitude <- getAmplitude(dsigma, fns, gs, points, ...)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / W^4) * abs(amplitude)^2
}
