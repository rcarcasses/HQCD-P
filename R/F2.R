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
  #flog.trace('Calling predict.F2')
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
getBestGs.F2 <- function(f2, fns, extended = FALSE) {
  f2Data <- expVal(f2)
  w <- 1/expErr(f2)^2
  # the -1 indicate that we don't want an intercept
  fit <- lm(f2Data ~ .-1, data = cbind(f2Data, fns), weights = w)
  #cat('getBestGs.F2', fit$coefficients, '\n')
  if(extended) {
    fit$gs <- fit$coefficients
    fit$value <- sum((fit$residuals^2) * w)
    fit
  }
  else
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
getFns.F2 <- function(f2, spectra, points) {
  alpha <- attr(f2, 'alpha')
  # precompute the modes to speed up
  mclapply(unique(points$Q2), getU1NNMode, alpha = alpha, mc.cores = cores)
  #flog.trace('Calling getFns.F2, alpha %s', alpha)
  spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
  reducer <- function(fNfun,  alpha = 0, extra = '')
    # iterate over each kernel's spectrum
    Reduce(function(acc, s) {
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      cbind(acc, Reduce(function(accspec, spec) {
        fn <- unlist(mclapply(apply(points, 1, as.list), function(row) {
          # we need to initialize the computation on each node
          init()
          fNfun(f2, row$Q2, row$x, spec$js, spec$wf, alpha)
        }, mc.cores = cores))
        val <- as.data.frame(cbind(accspec, fn))
        colnames(val)[length(val)] <- paste0(spec$name, extra)
        val
      }, s, init = c()))
  }, spectraForTZero, init = c(NA))

  df <- reducer(fN, alpha)
  if(alpha != 0)
    df <- cbind(df, reducer(fNNMC, alpha, '.NMC'))

  # remove the unneeded columns and return the result
  df[-which(names(df) == 'acc')]
}

#' @export
fN.F2 <- function(f2, Q2, x, J, wf, alpha) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getU1NNMode(Q2 = Q2, alpha = alpha)$factor # this is a spline fun
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  factor <- 1 / (tan(pi * J / 2) * 4 * pi^2 * 1 / 137)
  factor *  x^(1 - J) * Q2^J * integral$value
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
  factor <- 1 / (tan(pi * J / 2) * 4 * pi^2 * 1 / 137)
  factor * x^(1 - J) * Q2^J * integral$value
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

# this function attempts to reconstruct the wavefunctions from the raw data
reconstruct.F2 <- function(f2, js = c(1.3, 1.09)) {
  # get the data
  data <- cbind(list(F2 = expVal(f2)), expKinematics(f2), list(err = expErr(f2)))
  # get all the different values of Q2, once
  Q2s <- unique(data$Q2)
  Filter(function(q) !any(is.na(q)),
    lapply(Q2s, function(Q2) {
      # get the data of this Q2
      dataForQ2 <- data[data$Q2 == Q2,]
      # create the data which will be linearly fitted
      lmData <- cbind(list(F2 = dataForQ2$F2), as.data.frame(Reduce(function(acc, J) {
        # create the name of the column
        prop <- paste0('js', J)
        # assign the values
        acc[[prop]] <- apply(dataForQ2, 1, function(row) {
          x <- as.list(row)$x
          x^(1 - J)
        })
        # remove the names to prevent misbehaving
        names(acc[[prop]]) <- NULL
        acc
      }, js, init = list())))
      fit <- lm(F2 ~ . - 1, data = lmData, weights = dataForQ2$err^(-2))
      c(list(Q2 = Q2, rss = sum(fit$residuals^2)), fit$coefficients)
    }))
}

getBestExponents.F2 <- function(f2) {
  fOptim <- function(pars) {
    val <- sum(unlist(lapply(reconstruct.F2(f2, pars), `[[`, 'rss')))
    cat('pars', pars, ' value', val, '\n')
    val
  }
  optim(c(1.3, 1.1), fOptim)
}

#' For F2 the enlargement of the experimental data produces
#' some extra points with different values of x for each
#' different value of Q2
#' @export
enlargeData.F2 <- function(f2, ...) {
  # interesting values of x for producing a plot, for instace.
  xvals <- 10^seq(-6, -2, len = 50)
  Q2s <- sort(unique(f2$data$Q2))
  # for each one of the 34 previous values combine them with the xvals
  Reduce(rbind, lapply(Q2s, function(Q2) data.frame(Q2 = rep(Q2, length(xvals)), x = xvals)))
}

#' @export
plot.F2 <- function(f2, predicted) {
  data <- cbind(expKinematics(f2), list(F2 = expVal(f2)), list(err = expErr(f2)))
  # prepare the plot
  plot(10, 10,
       log = 'x',
       las = 1,
       xlim = c(1e-6, 1e-2),
       ylim = c(0.1, 1.8),
       xlab = 'x', ylab = expression(F[2]),
       xaxt = 'n', yaxt = 'n')
  hTicks <- seq(-6, -2, 1)
  vTicks <- seq(0, 2, 0.5)
  hLabels <- sapply(hTicks, function(i) as.expression(bquote(10^ .(i))))
  abline(v = 10^hTicks, h = vTicks, col = 'gray90', lty = 3)
  axis(1, at = 10^hTicks, labels = hLabels)
  axis(2, at = vTicks)
  numCols <- 4
  i <- 1
  cols <- sapply(seq(0.2, 0.7, len = numCols), gray)
  # get all different Q2s
  allQ2s <- sort(unique(f2$data$Q2))
  invisible(lapply(allQ2s, function(Q2) {
    # get the subset of data for the given Q2
    dataForQ2 <- data[data$Q2 == Q2,]
    colIndex <- mod(i, numCols) + 1
    # draw the experimental points
    lines(dataForQ2$x, dataForQ2$F2, type = 'p', pch = 19, col = cols[colIndex], cex = 0.7)
    # draw the experimental error bars
    with(dataForQ2,
      arrows(x, F2 - err, x, F2 + err, length = 0.02, angle = 90, code = 3, col = cols[colIndex])
    )
    # plot the predictions
    predictedForQ2 <- predicted[predicted$Q2 == Q2,]
    with(predictedForQ2, {
      lines(x, predicted, col = cols[colIndex])
    })
    # put the Q2 label
    textPosX <- 0.5 * min(dataForQ2$x)
    textPosY <- splinefun(predictedForQ2$x, predictedForQ2$predicted)(textPosX)
    boxed.labels(textPosX, textPosY, labels = paste(Q2), col = cols[colIndex],
                 cex = 0.8, xpad = 1.3, ypad = 1.3, border = T, bg = 'white')
    i <<- i + 1
  }))
}
