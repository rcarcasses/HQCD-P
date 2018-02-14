Sigma <- function(procName) {
  # add the generic DSigma class
  obs <- ProcessObservable(procName)
  class(obs) <- append(class(obs), 'Sigma', after = length(class(obs)) - 1)
  obs
}

# Here we need to take the t values from t = -1 to t = 0
#' @export
getNeededTVals.Sigma <- function(x) unique(expKinematics(x)$t)

#' Predicts the values of F2 for the points passed
#' @param sigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.Sigma <- function(sigma, fns, gs, points, ...) {
  tList <- seq(-1.0,0.0,0.05)
  gts <- apply(gs, 1, function(row) {
    g0 <- row[1]
    g1 <- row[2]
    g0 + g1 * tList
  })
  # compute the amplitude as a function of t
  amplitude <- rowSums(fns * gts, na.rm = TRUE)
  # get the Ws
  W <- points$W
  # return the differential cross section as a function of t
  dsigma <- (1 / W^4) * abs(amplitude)^2
  # We now need to compute the cross section as the integral over [-1,0] of dsigma
  dsigma <- splinfun(tlist, dsigma)
  integral <- integrate(function(x)  dsigma(x), -1.0, 0.0, stop.on.error = FALSE)
  integral$value
}

#' @export
getFns.Sigma <- function(sigma, points, spectra) {
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s[names(s) == ''],
                                           function(spec) paste0('fn.', spec$name)))))
  df <- data.frame(row.names = fnNames)
  as.data.frame(t(apply(points, 1, function(row) {
    W  <- row[1]
    Q2 <- row[2]
    t  <- seq(-1.0,0.0,0.05)
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    r <- unlist(lapply(spectraForT, function(s) {
      # s: spectrum of a single kernel, have many reggeons
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      lapply(s[names(s) == ''], function(spec) {
        fN(sigma, W, Q2, spec$js, spec$wf)
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  })))
}

fN.Sigma <- function(sigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalU1wf(sigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  (1 - 1i/ tan(pi * J / 2)) * W^(2*J) *integral$value
}

getExternalU1wf <- function(x, ...) UseMethod('getExternalU1wf')

getExternalU1wf.default <- function(x, ...) 'getExternalU1wf have to be implemented for this process'

#' @export
expVal.Sigma <- function(sigma) sigma$data$sigma
#' @export
expErr.Sigma <- function(sigma) sigma$data$deltaSigma
#' @export
expKinematics.Sigma <- function(sigma) sigma$data[c('Q2', 'W')]
