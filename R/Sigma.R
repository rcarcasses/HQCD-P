#' @export
Sigma <- function(procName) {
  # add the generic DSigma class
  obs <- DSigma(procName)
  class(obs) <- append(class(obs), 'Sigma', after = 2)
  obs
}

# Here we need to take the t values from t = -1 to t = 0
#' @export
getNeededTVals.Sigma <- function(x) seq(-1, 0, len = 20)

#' Predicts the values of F2 for the points passed
#' @param sigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.Sigma <- function(sig, fns, gs, points, ...) {
  ts <- getNeededTVals(sig)
  computeDSigma <- function(p) {
    fArgs <- list(...)
    # For each value of Q2 and W we need to insert many different values of t
    fArgs$points  <- cbind(Reduce(function(r1,i) rbind(r1, p), 1:(length(ts)), init = c()), ts)
    fArgs$fns     <- fns
    fArgs$gs      <- gs
    fArgs$generic <- 'predict'
    fArgs$object  <- sig
    do.call(NextMethod, fArgs)
  }
  dsigma <- apply(points, 1, computeDSigma)

  # We now need to compute the cross section as the integral over [-1,0] of dsigma
  dsigma <- splinfun(ts, dsigma)
  integral <- integrate(function(x)  dsigma(x), -1.0, 0.0, stop.on.error = FALSE)
  integral$value
}

#' @export
expVal.Sigma <- function(s) s$data$sigma
#' @export
expErr.Sigma <- function(s) s$data$deltaSigma
#' @export
expKinematics.Sigma <- function(s) s$data[c('Q2', 'W')]

