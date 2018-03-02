DSigma <- function(procName) {
  # add the generic DSigma class
  obs <- ProcessObservable(paste0(procName, 'DSigma'))
  class(obs) <- append(class(obs), 'DSigma', after = length(class(obs)) - 1)
  obs
}

#' @export
getNeededTVals.DSigma <- function(x) unique(expKinematics(x)$t)

getAmplitude <- function(x, ...) UseMethod('getAmplitude')
getAmplitude.DSigma <- function(dsigma, fns, gs, points, ...) {
  # compute g(t) for each corresponding fn, this return a dataframe
  # were each column is the value (a vector) of g(t) for the given values
  # of t and the n column is the g(t) of the n fn in the fns dataframe
  # see the test file for some explanation 'tests/testthat/test_DSigma.R'
  gts <- apply(gs, 1, function(row) {
    rowSums(t(row * t(outer(t, 0:(length(gs) - 1), `^`))))
  })
  rowSums(fns * gts, na.rm = TRUE)
}

#' Get fns times dJdt
#' @export
getFns.DSigma <- function(dsigma, spectra, points) {
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s, function(spec) paste0('fn.', spec$name)))))
  df <- data.frame(row.names = fnNames)
  as.data.frame(t(apply(points, 1, function(row) {
    row <- as.list(row)
    Q2 <- row$Q2
    W  <- row$W
    t  <- row$t
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    r <- unlist(lapply(spectraForT, function(s) {
      # s: spectrum of a single kernel, have many reggeons
      # iterate over each Reggeon for the given spectrum
      lapply(s, function(spec) {
        spec$dJdt * fN(dsigma, W, Q2, spec$js, spec$wf)
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  })))
}

fN.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalStateFactor(dsigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  (1 - 1i/ tan(pi * J / 2)) * W^(2*J) *integral$value
}

getExternalStateFactor <- function(x, ...) UseMethod('getExternalStateFactor')

getExternalStateFactor.default <- function(x, ...) 'getExternalStateFactor have to be implemented for this process'

#' @export
expVal.DSigma <- function(dsigma) dsigma$data$dsigma
#' @export
expErr.DSigma <- function(dsigma) dsigma$data$deltaDSigma
#' @export
expKinematics.DSigma <- function(dsigma) dsigma$data[c('Q2', 'W', 't')]
