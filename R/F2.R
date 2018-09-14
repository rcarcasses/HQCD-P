#' @export
F2 <- function() ProcessObservable('F2')

#' For DIS we only need t = 0 value
#' @export
getNeededTVals.F2 <- function(f2) c(0)

#' Predicts the values of F2 for the points passed
#' @param f2 the object over which the prediction will happend
#' @param Izs the return of the getIzs function on this object. Contains
#' the values of the integrals for each one of the Reggeons involved
#' @param IzsBar the return of the getIzBars function on this object
#' the observable values.
#' @param ... Additional arguments will be ignored.
#' @export
predict.F2 <- function(f2, Izs, IzsBar, ...) {
  # get the predicted vector: sum_i(I_i * Ibar_i)
  val <- rowSums(Izs * IzsBar, na.rm = TRUE)
  val
}

#' This function returns a term in the sum for F2 with the
#' exception of the overall constant.
#' @param points The kinematical points over which the Izs will be computed
#' @param spectra An object with the spectrum of all the kernels for different
#' values of t required for computing the integrals
#' @return A data frame where the first columns are the points passed
#' and the next ones are the values of Izs for each one of the kernels
#' @export
getIzs.F2 <- function(f2, spectra, points) {
  alpha <- attr(f2, 'alpha')
  # precompute the modes to speed up
  mclapply(unique(points$Q2), getU1NNMode, alpha = alpha, mc.cores = cores)
  #flog.trace('Calling getIzs.F2, alpha %s', alpha)
  spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
  reducer <- function(IzNfun,  alpha = 0, extra = '')
    # iterate over each kernel's spectrum
    Reduce(function(acc, s) {
      # iterate over each Reggeon for the given spectrum
      # remember, the tr1 and tr2 are not data about the reggeons
      cbind(acc, Reduce(function(accspec, spec) {
        fn <- unlist(mclapply(apply(points, 1, as.list), function(row) {
          # we need to initialize the computation on each node
          init()
          IzNfun(f2, row, spec, alpha)
        }, mc.cores = cores))
        # stop in case of an error, and display some information
        if(is.character(fn[1])) {
          flog.warn('Non numeric result in F2')
          flog.warn('%s ', fn[1])
          stop()
        }
        val <- as.data.frame(cbind(accspec, fn))
        colnames(val)[length(val)] <- paste0(spec$name, extra)
        val
      }, s, init = c()))
  }, spectraForTZero, init = c(NA))

  df <- if('IzN' %in% names(attributes(f2)))
    reducer(attr(f2, 'IzN'), 0)
  else
    reducer(IzN, 0)
  # if considering non minimal coupling
  if(alpha != 0)
    df <- cbind(df, reducer(IzNNMC, alpha, '.NMC'))
  # remove the unneeded columns and return the result
  df[-which(names(df) == 'acc')]
}

#' @export
IzN.F2 <- function(f2, kin, spec, alpha) {
  Q2 <- kin$Q2
  x  <- kin$x
  J  <- spec$js
  wf <- spec$wf
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getU1NNMode(Q2 = Q2, alpha = alpha)$factor # this is a spline fun
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  factor <- 1 / (4 * pi^2 * 1 / 137)
  factor *  x^(1 - J) * Q2^J * integral$value
}

IzNNMC.F2 <- function(f2, kin, spec, alpha) {
  Q2 <- kin$Q2
  x  <- kin$x
  J  <- spec$js
  wf <- spec$wf
  mode  <- getU1NNMode(Q2 = Q2, alpha)
  fQ2   <- mode$fQ2
  dfQ2  <- mode$dfQ2
  # here is where the non-minimal coupling changes are implemented
  # we have included the factor of exp(-2 As) in the definition of the Ds
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  # contribution from the transverse part
  t2fun   <- splinefun(z, fQ2(z) * DperpPsi(wf))
  # contribution from the longitudinal part
  t3fun   <- splinefun(z, (dfQ2(z) / Q2) * DparallelPsi(wf))
  integral <- integrate(function(x) t1fun(x) * (t2fun(x) + t3fun(x)), z[1], z[length(z)], stop.on.error = FALSE)
  factor <- 1 / (4 * pi^2 * 1 / 137)
  factor * x^(1 - J) * Q2^J * integral$value
}

#' @export
getIzsBar.F2 <- function(f2, spectra, points, zstar, hpars) {
  alpha <- attr(f2, 'alpha')
  #flog.trace('Calling getIzs.F2, alpha %s', alpha)
  spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
  reducer <- function(IzNBarfun,  alpha = 0, extra = '')
    # iterate over each kernel's spectrum
    Reduce(function(acc, s) {
      # iterate over each Reggeon for the given spectrum
      cbind(acc, Reduce(function(accspec, spec) {
        # repeat the results as needed
        fn <- rep(
                Im(IzNBarfun(f2, NULL, spec, zstar, hpars)),
              length(points[[1]]))
        val <- as.data.frame(cbind(accspec, fn))
        colnames(val)[length(val)] <- paste0(spec$name, extra)
        val
      }, s, init = c()))
    }, spectraForTZero, init = c(NA))
  # use the attribute if defined while computing the IzNBar integrals
  df <- if('IzNBar' %in% names(attributes(f2)))
    reducer(attr(f2, 'IzNBar'), 0)
  else
    reducer(IzNBar, 0)
  # TODO: implement the IzNNMCBar
  #if(alpha != 0)
  #  df <- cbind(df, reducer(IzNNMCBar, alpha, '.NMC'))

  # remove the unneeded columns and return the result
  df[-which(names(df) == 'acc')]
}

#' compute Dperp acting on psi
#' @export
DperpPsi <- function(wf) {
  wffun <- splinefun(wf$x, wf$y)
  wfder1 <- wffun(z, deriv = 1)
  # return object
  exp(-2 * As) * (
    Asder1 * wfder1 +
    (0.5 * Asder1^2 + Asder1 * Phider1) * wffun(z)
  )
}
#' compute Dparallel acting on psi
#' @export
DparallelPsi <- function(wf) {
  wffun <- splinefun(wf$x, wf$y)
  wfder1 <- wffun(z, deriv = 1)
  wfder2 <- wffun(z, deriv = 2)
  # return object
  exp(-2 * As) * (
    wfder2
    + 2 * (Phider1 - Asder1) * wfder1
    + (0.75 * Asder1^2 - 2 * Asder1 * Phider1 + Phider1^2 + 0.5 * Asder2 + Phider2) * wffun(z)
  )
}


#' @export
expVal.F2 <- function(f2) f2$data$F2
#' @export
expErr.F2 <- function(f2) f2$data$err
#' @export
expKinematics.F2 <- function(f2) as.data.frame(f2$data[c('Q2', 'x')])

# this function attempts to reconstruct the wavefunctions from the raw data
#' @export
reconstruct.F2 <- function(f2, js = c(1.17, 1.09), getPsiNs = FALSE, isAdS = FALSE) {
  # get the data
  data <- cbind(list(F2 = expVal(f2)), expKinematics(f2), list(err = expErr(f2)))
  # get all the different values of Q2, once
  Q2s <- unique(data$Q2)
  df <- data.frame(Q2 = c(), rss = c(), c1 = c(), c2 = c())
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
      # get the standard errors
      s.e <- summary(fit)$coefficients[,2]
      df <<- rbind(df, data.frame(Q2 = Q2,
                                  rss = sum(fit$residuals^2),
                                  c1 = fit$coefficients[1],
                                  dc1 = s.e[1],
                                  c2 = fit$coefficients[2],
                                  dc2 = s.e[2]))
    }))
  rownames(df) <- c()
  fNs <- na.omit(df)
  r <- if(getPsiNs)
    list(fNs = fNs, psiNs = reconstructPsiNs(fNs, js, isAdS))
  else
    fNs
  # mark this as an object that contains reconstruction data of F2
  class(r) <- append('F2reconstructed', class(r))
  # pass the used js as attribute
  attr(r, 'js') <- js
  r
}

#' @export
plot.F2reconstructed <- function(r, method = 'fNs', ...) {
  js <- attr(r, 'js')
  print(js)
  df <- if(is.null(r$fNs)) r else r$fNs
  f0 <- function(pch = 20, col = 'black') {
    # plot f0
    lines(log(df$Q2), df$c1, pch = pch, col = col, type = 'p')
    # plot the error bars with 2 sigma
    arrows(log(df$Q2), df$c1 - 2 * df$dc1, log(df$Q2), df$c1 + 2 * df$dc1,
           length = 0.04, angle = 90, code = 3, col = col)
    # plot a background grid to improve readability
    abline(h = 0.5 * (0:3), v = 2 * (-1:3), lty = 3, col = 'lightgrey')
    # plot an approximating function
    f0.lm <- lm(c1 ~ I(log(Q2))+I(log(Q2)^2), data = df)
    logQ2s <- seq(-2, 6, len = 200)
    lines(logQ2s, predict(f0.lm, newdata = data.frame(Q2 = exp(logQ2s))), col = col)
  }
  f1 <- function(pch = 20, col = 'black') {
    # plot f1
    lines(log(df$Q2), df$c2, pch = pch, col = col, type = 'p')
    # plot the error bars with 2 sigma
    arrows(log(df$Q2), df$c2 - 2 * df$dc2, log(df$Q2), df$c2 + 2 * df$dc2,
           length = 0.04, angle = 90, code = 3, col = col)
    # plot a background grid to improve readability
    abline(h = 0.5 * (-4:1), v = 2 * (-1:3), lty = 3, col = 'lightgrey')
    # plot an approximating function
    f1.lm <- lm(c2 ~ I(log(Q2))+I(log(Q2)^2), data = df)
    logQ2s <- seq(-2, 6, len = 200)
    lines(logQ2s, predict(f1.lm, newdata = data.frame(Q2 = exp(logQ2s))), col = col)
  }

  switch (method,
    'fNs' = {
      plot(-100, xlim = c(-2, 6), ylim = c(-2, 2),
           ylab = expression(f[0] ~ ', '~ f[1]),
           xlab = expression(log(Q^2)),
           main = eval(substitute(expression(f[0](Q^2) ~ ' and ' ~ f[1](Q^2) ~ ' for ' ~ j[0] == j0 ~ ' and ' ~ j[1] ==  j1), list(j0 = format(js[1], digits = 4), j1 = format(js[2], digits = 4))))
           )
      f0(col = 'red')
      f1(col = 'blue')
    },
    'f0'  = {
      plot(-100, xlim = c(-2, 6), ylim = c(-0.3, 2),
           ylab = expression(f[0](Q^2)),
           xlab = expression(log(Q^2)),
           main = expression(f[0](Q^2)))
      f0()
    },
    'f1'  = {
      plot(1, xlim = c(-2, 6), ylim = c(-2, 0.3),
           ylab = expression(f[1](Q^2)),
           xlab = expression(log(Q^2)),
           main = expression(f[1](Q^2)))
      f1()
    }
  )
}



# From the tests I have made it seems like the most important
# ingredient to get the right shapes for the kernel wavefuntions
# is the dilaton profile: you can run AdS and warped AdS cases
# with zero dilaton and the reconstruction does not work well,
# however, if you consider a non trivial dilaton in both cases
# one get reasonable shapes, particularly better for the warped
# AdS case.
reconstructPsiNs <- function(fNs, js = c(1.17, 1.09), isAdS = FALSE) {
  if(isAdS) {
    Asfun <- splinefun(z, -log(z))
    Phifun <- splinefun(z, rep(0, length(z)))
  } else {
    Asfun <- splinefun(z, As)
    Phifun <- splinefun(z, Phi)
  }
  fNs <- cbind(fNs, z = 1/sqrt(r$Q2))
  lapply(1:length(js), function(jsIndex) {
    as.data.frame(cbind(z = 1/sqrt(r$Q2),
          psi = apply(r, 1, function(row) {
            row <- as.list(row)
            exp((js[jsIndex] - (1/2)) * Asfun(1/row$z)) * exp(-Phifun(row$z)) * row$Q2^(-js[jsIndex]) * row[[paste0('c', jsIndex)]]
          })))
  })
}

#' @export
plotReconstructedf0f1 <- function(js = c(1.17, 1.09), isAdS = FALSE) {
  r <- reconstruct(F2(), js)
  plot(r$Q2, r$c1, log = 'x', pch = 20, col = 'blue',
       ylim = c(-1.5, 1.5),
       ylab = expression(paste(f[0], ', ', f[1])),
       xlab = expression(Q^2))
  lines(r$Q2, r$c2, type = 'p', pch = 20, col = 'red')
  abline(h = 0.5 * (-3:3), v = c(0.2, 2, 20, 200), col = 'gray80', lty = 3)
}

#' @export
plotReconstructedPsiNs <- function(js = c(1.17, 1.09), isAdS = FALSE) {
  psin <- reconstructPsiNs(js, isAdS)
  plot(psin[[1]]$z, psin[[1]]$psi, col = 'red',
       ylim = c(-0.25, 0.5),
       ylab = expression(paste(psi[0], ', ', psi[1])),
       xlab = expression(z), pch = 20)
  lines(psin[[2]]$z, psin[[2]]$psi, col = 'blue', type = 'p', pch = 20)
  abline(h = 0.1 * (-3:5), v = 0.5 * (0:5), col = 'gray80', lty = 3)
}

#' @export
getBestExponents.F2 <- function(f2) {
  fOptim <- function(pars) {
    val <- sum(reconstruct(f2, c(pars, 1.09))$rss)
    cat('pars', pars, ' value', val, '\n')
    val
  }
  optim(c(1.3), fOptim)
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
       xlim = c(9e-7, 1e-2),
       ylim = c(0.1, 1.8),
       xlab = 'x', ylab = expression(F[2]),
       xaxt = 'n', yaxt = 'n')
  hTicks <- seq(-6, -2, 1)
  vTicks <- seq(0, 2, 0.5)
  hLabels <- sapply(hTicks, function(i) as.expression(bquote(10^ .(i))))
  abline(v = 10^hTicks, h = vTicks, col = 'gray90', lty = 3)
  axis(1, at = 10^hTicks, labels = hLabels)
  axis(2, at = vTicks, las = 1)
  numCols <- 3
  i <- 1
  cols <- c('blue', 'red', 'green')#sapply(seq(0.2, 0.7, len = numCols), gray)
  # get all different Q2s
  allQ2s <- sort(unique(f2$data$Q2))
  invisible(lapply(allQ2s, function(Q2) {
    # get the subset of data for the given Q2
    dataForQ2 <- data[data$Q2 == Q2,]
    colIndex <- mod(i, numCols) + 1
    # draw the experimental points
    lines(dataForQ2$x, dataForQ2$F2, type = 'p', pch = 1, col = cols[colIndex], cex = 0.7)
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
    textPosX <- 1.1 * textPosX
    boxed.labels(textPosX, textPosY, labels = paste(Q2), col = cols[colIndex],
                 cex = 0.8, xpad = 1.3, ypad = 1.3, border = FALSE, bg = 'white', xlog = TRUE)
    i <<- i + 1
  }))
}
