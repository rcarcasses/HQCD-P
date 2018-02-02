#' @export
DVCSDSigma <- function() ProcessObservable('DVCSDSigma')

#' @export
getNeededTVals.DVCSDSigma <- function(x) unique(expKinematics(x)$t)

#' Predicts the values of F2 for the points passed
#' @param dvcsDSigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt for dvcs will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.DVCSDSigma <- function(dvcsDSigma, points, gs, spectra) {
  # compute the amplitude
  amplitude <- computeAmplitude(dvcsDSigma, points, gs, spectra)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / W^4) * abs(amplitude)^2
}

#' Compute the amplitude
#' @return A data frame where the first columns are the points passed
#' and the next ones are the values of fns for each one of the kernels
#' @export
computeAmplitude.DVCSDSigma <- function(dvcsDSigma, points, gs, spectra) {
  apply(points, 1, function(row) {
    W  <- row[1]
    Q2 <- row[2]
    t  <- row[3]
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    sum(unlist(mapply(function(s, gker) {
      # s: spectrum of a single kernel, have many reggeons
      # gker: list of coefficients for the couplings for a given kernel
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      sum(unlist(mapply(function(spec, g) {
          # spec: spectrum of a given Reggeon inside this kernel
          # g: the associated coupling parameters
          # we are approximating g(t) = g0 + g1 * t + ...
          (g$g0  + g$g1 * t) * fN(dvcsDSigma, W, Q2, spec$js, spec$wf)
        }, s[names(s) == ''], gker)))
      }, spectraForT, gs)))
  })
}


fN.DVCSDSigma <- function(dvcsDSigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getU1NNMode(Q2 = Q2)$fQ
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  (1 - 1i/ tan(pi * J / 2)) * W^(2*J) *integral$value
}

#' @export
expVal.DVCSDSigma <- function(dvcs) dvcs$data$dsigma
#' @export
expErr.DVCSDSigma <- function(dvcs) dvcs$data$deltaDSigma
#' @export
expKinematics.DVCSDSigma <- function(dvcs) dvcs$data[c('Q2', 'W', 't')]
