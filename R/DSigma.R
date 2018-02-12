DSigma <- function(procName) {
  # add the generic DSigma class
  obs <- ProcessObservable(procName)
  class(obs) <- append(class(obs), 'DSigma', after = length(class(obs)) - 1)
  obs
}

#' @export
getNeededTVals.DSigma <- function(x) unique(expKinematics(x)$t)

#' Predicts the values of F2 for the points passed
#' @param dsigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.DSigma <- function(dsigma, fns, gs, points, ...) {
  gts <- apply(gs, 1, function(row) {
    g0 <- row[1]
    g1 <- row[2]
    g0 + g1 * points$t
  })
  # compute the amplitude
  amplitude <- rowSums(fns * gts, na.rm = TRUE)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / W^4) * abs(amplitude)^2
}

#' @export
getFns.DSigma <- function(dsigma, points, spectra) {
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s[names(s) == ''],
                                           function(spec) paste0('fn.', spec$name)))))
  df <- data.frame(row.names = fnNames)
  as.data.frame(t(apply(points, 1, function(row) {
    W  <- row[1]
    Q2 <- row[2]
    t  <- row[3]
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    r <- unlist(lapply(spectraForT, function(s) {
      # s: spectrum of a single kernel, have many reggeons
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      lapply(s[names(s) == ''], function(spec) {
        fN(dsigma, W, Q2, spec$js, spec$wf)
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  })))
}

fN.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalU1wf(dsigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  (1 - 1i/ tan(pi * J / 2)) * W^(2*J) *integral$value
}

getExternalU1wf <- function(x, ...) UseMethod('getExternalU1wf')

getExternalU1wf.default <- function(x, ...) 'getExternalU1wf have to be implemented for this process'

gradG.DSigma <- function(dsigma, fns, gs) {
  points <- expKinematics(dsigma)
  gts <- apply(gs, 1, function(row) {
    g0 <- row[1]
    g1 <- row[2]
    g0 + g1 * points$t
  })
  # keep in mind that * is not matrix multiplication
  df  <- as.data.frame(2 * Re (as.matrix(fns) * rowSums(Conj(as.matrix(fns)) * gts)))
  dft <- df * points$t
  cbind(df, dft)
}

#' @export
expVal.DSigma <- function(dsigma) dsigma$data$dsigma
#' @export
expErr.DSigma <- function(dsigma) dsigma$data$deltaDSigma
#' @export
expKinematics.DSigma <- function(dsigma) dsigma$data[c('Q2', 'W', 't')]