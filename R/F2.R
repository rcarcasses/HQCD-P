#' @export
F2 <- function() ProcessObservable('F2')

#' For DIS we only need t = 0 value
#' @export
getNeededTVals.F2 <- function(f2) c(0)

#' Predicts the values of F2 for the points passed
#' @param f2 the object over which the prediction will happend
#' @param fns the return of the getFns function on this object. Contains
#' the values of the integrals for each one of the Reggeons involved
#' @param gs the coefficients, which together with the fns, can be used to compute
#' the observable values.
#' @param ... Additional arguments will be ignored.
#' @export
predict.F2 <- function(f2, fns, gs, ...) {
  # if we are passing the full gs then use only the g(t=0) part
  if(is.data.frame(gs))
    gs <- gs[,1]

  # get the predicted vector: sum_i(g_i * fn_i)
  val <- rowSums(mapply(function(c, g) c * g, fns, gs), na.rm = TRUE)

  if(is.null(attr(f2, 'complete')))
    val
  else
    list(val = val, coefficients = gs)
}

#' Given a pre-computed set of fns find the best coefficients
#' such that prediction matches best with the F2 data.
#' In this case the dependence of F2 with the gs is linear, which
#' makes the problem very simple.
#' @export
getBestGs.F2 <- function(f2, fns) {
  f2Data <- expVal(f2)
  w <- 1/expErr(f2)^2
  # the -1 indicate that we don't want an intercept
  fit <- lm(f2Data ~ .-1, data = cbind(f2Data, fns), weights = w)
  #cat('getBestGs.F2', fit$coefficients, '\n')
  fit$coefficients
}

#' This function returns a term in the sum for F2 with the
#' exception of the overall constant.
#' @param points The kinematical points over which the fns will be computed
#' @param spectra An object with the spectrum of all the kernels for different
#' values of t required for computing the integrals
#' @return A data frame where the first columns are the points passed
#' and the next ones are the values of fns for each one of the kernels
#' @export
getFns.F2 <- function(f2, points, spectra) {
  spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
  reducer <- function(alpha = 0)
    # iterate over each kernel's spectrum
    Reduce(function(acc, s) {
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      val0  <- cbind(acc, Reduce(function(accspec, spec) {
        fn <- apply(points, 1, function(row) {
          Q2 <- row[1]
          x  <- row[2]
          if(alpha == 0)
            fN(f2, Q2, x, spec$js, spec$wf)
          else
            fNNMC(f2, Q2, x, spec$js, spec$wf, alpha)
        })
        val <- as.data.frame(cbind(accspec, fn))
        colnames(val)[length(val)] <- if(alpha == 0) spec$name else paste0(spec$name, '.NMC')
        val
      }, s, init = c()))
      val0
  }, spectraForTZero, init = c(NA))

  df <- reducer()
  if(attr(f2, 'alpha') != 0)
    df <- cbind(df, reducer(attr(f2, 'alpha')))

  # remove the unneeded columns and return the result
  df[-which(names(df) == 'acc')]
}

#' @export
fN.F2 <- function(f2, Q2, x, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getU1NNMode(Q2 = Q2, alpha = 0)$factor # this is a spline fun
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  x^(1 - J) * Q2^J * integral$value
}

fNNMC.F2 <- function(f2, Q2, x, J, wf, alpha) {
  mode  <- getU1NNMode(Q2 = Q2, alpha)
  fQ2   <- mode$fQ2
  dfQ2  <- mode$dfQ2
  # here is where the non-minimal coupling changes are implemented
  # we have included the factor of exp(-2 As) in the definition of the Ds
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  # contribution from the transverse part
  t2fun   <- splinefun(z, fQ2(z) * DperpPsi(wf))
  # contribution from the longitudinal part
  t3fun   <- splinefun(z, (dfQ2(z) / Q2) * DparallelPsi(wf, J))
  integral <- integrate(function(x) t1fun(x) * (t2fun(x) + t3fun(x)), z[1], z[length(z)], stop.on.error = FALSE)
  # cat('int ', integral$value, '\n')
  x^(1 - J) * Q2^J * integral$value
}

#' compute Dperp acting on psi
#' @export
DperpPsi <- function(wf) {
  wffun <- splinefun(wf$x, wf$y)
  wfder1 <- wffun(z, deriv = 1)
  # return object
  exp(-2 * As) * (
    (Asder1 * Phider1 - 1.5 * Asder1^2) * wffun(z) + Asder1 * wfder1
  )
}
#' compute Dparallel acting on psi
#' @export
DparallelPsi <- function(wf, J) {
  wffun <- splinefun(wf$x, wf$y)
  wfder1 <- wffun(z, deriv = 1)
  wfder2 <- wffun(z, deriv = 2)
  # return object
  exp(-2 * As) * (
    wfder2 + 2 * (Phider1 - Asder1) * wfder1 +
     (Phider2 + (J - 2.5) * Asder2 - Asder1 * (2 * Phider1 + J - 1)
      + Phider1^2 + 0.75 * Asder1^2) * wffun(z)
    )
}

# DEPRECATED
gradG.F2 <- function(obs, fns, gs) {
  # derivatives respect of g1 are zero
  cbind(fns, as.data.frame(matrix(0, ncol = length(fns), nrow = length(fns[[1]]))))
}

#' @export
expVal.F2 <- function(f2) f2$data$F2
#' @export
expErr.F2 <- function(f2) f2$data$err
#' @export
expKinematics.F2 <- function(f2) as.data.frame(f2$data[c('Q2', 'x')])
