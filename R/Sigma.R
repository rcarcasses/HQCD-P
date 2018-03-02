#' @export
Sigma <- function(procName, tmin, tmax) {
  # add the generic DSigma class
  obs <- ProcessObservable(paste0(procName, 'Sigma'))
  class(obs) <- append(class(obs), 'Sigma', after = 2)
  attr(obs, 'tmin') <- tmin
  attr(obs, 'tmax') <- tmax
  # create an internal dsigma object and
  # enhance artificially the data including N values of t for
  # each entry point. These will be required for the computation
  # of the integral
  obs$dsigma <- DSigma(procName)
  obs
}

#' @export
enhanceDataWithTs <- function(sigma, data) {
  ts <- getNeededTVals(sigma)
  as.data.frame(Reduce(function(acc, r) {
    val <- rbind(acc, r)
    rownames(val) <- NULL
    val
  }, unlist(
    apply(data, 1,
          function(row)
            lapply(ts, function(t) c(list(t = t), row))),
    recursive = FALSE)
  ), row.names = NULL)
}

# Here we need to take the t values from t = -1 to t = 0
#' @export
getNeededTVals.Sigma <- function(x) seq(attr(x, 'tmin'), attr(x, 'tmax'), len = 10)

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
  # compute all the needed differential cross-sections
  # For each value of Q2 and W we need to insert many different values of t
  sig$dsigma$data <- enhanceDataWithTs(sig, points)
  dsigma <- predict(sig$dsigma, fns = fns, gs = gs, ...)

  ts <- getNeededTVals(sig)
  blockSize <- length(points) / length(ts)
  sigma <- unlist(lapply(1:blockSize, function(i) {
    # get the dsigma data computed for the value of Q2 and W, etc...
    # in the present block
    ds <- dsigma[1:length(ts) + (i - 1) * length(ts)]

    # we now need to compute the cross section as the integral over [-1,0] of dsigma
    dsFun <- splinfun(ts, dsigma)
    integral <- integrate(function(x)  dsFun(x), attr(x, 'tmin'), attr(x, 'tmax'), stop.on.error = FALSE)
    integral$value
  }))
  sigma
}

#' @export
getFns.Sigma <- function(sigma, points, spectra) {
  # return the fns for the enhanced points of the correspondent dsigma object
  getFns(sigma$dsigma, spectra = spectra, points = enhanceDataWithTs(sigma, points))
}

#' @export
expVal.Sigma <- function(s) s$data$sigma
#' @export
expErr.Sigma <- function(s) s$data$deltaSigma
#' @export
expKinematics.Sigma <- function(s) s$data[c('Q2', 'W')]

