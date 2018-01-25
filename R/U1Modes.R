
#' @export
modesEnv <- new.env()

#' @export
getMode <- function(Q2, A0 = NULL, h, alpha = -0.06718486) {
  if(is.null(A0)) {
    cat('[WARN] A0 is null in getMode(), setting to the IHQCD value\n')
    A0 <- get('A0', envir = ihqcdEnv)
  }

  key <- paste('u1m', A0, h, alpha, sep = '-')
  if(!exists(key, envir = modesEnv)){  # compute them if needed
    if(!rredis::redisExists(key))
      findU1modesForAllQ2(key, A0, h, alpha)

    # copy all from redis to a local environment
    assign(key, rredis::redisGet(key), envir = modesEnv)
    # here we also compute all the factors
    # cat('Computing mode factors...')
    modes <- as.list(get(key, envir = modesEnv))
    factors <- list()
    for (i in 1:length(modes$Q2)) {
      m       <- list(x = modes$x[[i]], y = modes$y[[i]])
      f       <- m$y
      fspline <- splinefun(m$x, f)
      df      <- fspline(m$x, deriv = 1)
      t1      <- f^2 + (1/modes$Q2[[i]]) * df^2
      factor  <- splinefun(m$x, t1)
      # for the non minimal coupling case we need the following pieces
      f2Q     <- splinefun(m$x, f^2)
      f2Qder1 <- splinefun(m$x, df^2)
      Q2str   <- paste0(modes$Q2[[i]])
      factors[[Q2str]] <- list(factor = factor, f2Q = f2Q, f2Qder1 = f2Qder1)
    }
    assign('u1mfac', factors, envir = modesEnv)
    # cat('done.\n')
    # now update the local modes
    assign(key, modes, envir = modesEnv)
  }

  modes   <- as.list(get(key, envir = modesEnv))
  factors <- get('u1mfac', envir = modesEnv)

  for (i in 1:length(modes$Q2)) {
    if(modes$Q2[[i]] == Q2) {
      el      <- factors[[paste0(Q2)]]
      factor  <- el$factor
      f2Q     <- el$f2Q
      f2Qder1 <- el$f2Qder1

      return(list(x = modes$x[[i]],
                  y = modes$y[[i]],
                  factor = factor,
                  f2Q = f2Q,
                  f2Qder1 = f2Qder1
                ))
    }
  }
  cat(paste('[WARN] No mode found for Q^2=', Q2,'\n'))
}

# these are the non-normalizable modes
#' @export
findU1Modes <- function(Q2 = 3.5, alpha = 0) {
  flog.debug(paste('finding mode', Q2, 'alpha', alpha))
  x             <- get('z', envir = ihqcdEnv)
  As            <- get('As', envir = ihqcdEnv)
  Asder1        <- get('Asder1', envir = ihqcdEnv)
  Asder2        <- get('Asder2', envir = ihqcdEnv)
  Asder3        <- get('Asder3', envir = ihqcdEnv)
  Phider1       <- get('Phider1', envir = ihqcdEnv)
  Asfun         <- splinefun(x, As)
  Asder1fun     <- splinefun(x, Asder1)
  Asder2fun     <- splinefun(x, Asder2)
  Asder3fun     <- splinefun(x, Asder3)
  Phider1fun    <- splinefun(x, Phider1)
  # this is useful for the non minimal coupling case
  fact0fun      <- splinefun(x, (1 - 2 * exp(-2 * As) * alpha * Asder2) / (1 - 2 * exp(-2 * As) * alpha * Asder1^2))
  fact1fun      <- splinefun(x, exp(2 * As) - 2 * alpha * Asder2)

  # now we need to define the differential equation for the U(1) field modes
  fun <- function(x, y, pars) {
    dy2 <- Q2 * y[1] - (Asder1fun(x) - Phider1fun(x)) * y[2]
    dy1 <- y[2]
    list(c(dy1, dy2))
  }
  funAlpha <- function(x, y, pars) {
    # remark: we are using Q2=q^2 and not the usual -q^2 since we are using the GR metric convention for eta
    dy2 <- Q2 * y[1] * fact0fun(x) - (Asder1fun(x) - Phider1fun(x)
                        + 2 * alpha * (2 * Asder1fun(x) * Asder2fun(x) - Asder3fun(x)) / fact1fun(x)) * y[2]
    dy1 <- y[2]
    list(c(dy1, dy2))
  }
  odeFun <- fun
  if (abs(alpha) > 0)
    odeFun <- funAlpha
  # return the solution as a data frame
  as.data.frame(bvpcol(yini = c(y = 1, dy = NA), x = x, func = odeFun, yend = c(0, NA), nmax = 100000, atol = 1e-9))
}

#' @export
getU1NormalizableModesPotential <- function() {
  x          <- get('z', envir = ihqcdEnv)
  A          <- get('A', envir = ihqcdEnv)
  Ader1      <- get('Ader1', envir = ihqcdEnv)
  Ader2      <- get('Ader2', envir = ihqcdEnv)
  lambda     <- get('lambda', envir = ihqcdEnv)
  lambdader1 <- get('lambdader1', envir = ihqcdEnv)
  lambdader2 <- get('lambdader2', envir = ihqcdEnv)

  Phider1 <- lambdader1 / lambda
  Phider2 <- Phider1^2 - lambdader2 / lambda
  Atfun   <- splinefun(x, A - log(lambda) / 3)
  Atder1  <- Ader1 - Phider1 / 3
  Atder2  <- Ader2 - Phider2 / 3

  y <- 0.25 * Atder1^2 + 0.5 * Atder2
  list(x = x, y = y, Atfun = Atfun)
}

#' @export
findU1NormalizableModes <- function(n = 10) {
  p <- getU1NormalizableModesPotential()

  setPotential(p$x, p$y)
  computeSpectrum(n)
  wfs <- getWavefunctions()
  Qs <- getEnergies()
  mapply(function(wf, Q2) {
    list(x = wf$x, y = (abs(wf$y[100]) / wf$y[100]) * wf$y * exp(0.5 * p$Atfun(wf$x)), Q2 = Q2)
  }, wfs, Qs, SIMPLIFY = FALSE)
}

#' @export
plotU1NormalizableMass2 <- function(n = 20) {
  modes <- findU1NormalizableModes(n)
  Q2s <- sapply(modes, function(m) m$Q2)
  x <- 1:n
  plot(x, y = Q2s, ylab = expression('m'^2), xlab = 'n', xlim = c(0, n), ylim = c(0, max(Q2s)),
       pch=ifelse(x==1, 19, 1))
  label <- getParExpression('m1', format(sqrt(Q2s[1]), digits = 2))
  boxed.labels(1.8, Q2s[1], labels = label, cex = 1, xpad = 1, ypad = 1, border = FALSE, bg = 'transparent')

  return()
}

#' @export
plotU1NormalizableModes <- function(indices = c(1,3,10)) {
  n <- max(indices)
  scale <- 10
  modes <- findU1NormalizableModes(n)
  p <- getU1NormalizableModesPotential()

  plot.new()
  plot(p$x, p$y, ylim = c(-6, 12), xlim = c(0, 5), lwd = 3, col = 'gray', type = 'l', ylab = 'U', xlab = 'z')
  abline(h = c(0), v = c(0), col = 'gray', lwd = 2)

  lapply(indices, function(index) {
    # get the mode
    m <- modes[[index]]
    # get the position of the first maximum
    mfun <- splinefun(m$x, m$y)
    mfunder1 <- splinefun(m$x, mfun(m$x, deriv = 1))
    maxX <- uniroot.all(f = mfunder1, c(0, 2))[[1]]
    maxY <- scale * mfun(maxX)
    # plot it
    lines(m$x, scale * m$y)
    # and plot a small label indicating the mode index
    label <- getParExpression(paste0('m', index), format(sqrt(m$Q2), digits = 2))
    boxed.labels(maxX + 0.4, maxY + 0.5, labels = label, cex = 1, xpad = 1, ypad = 1, border = FALSE, bg = 'transparent')
  })
  return()
}

Q2s <- c(0.10, 0.15, 0.20, 0.25, 0.35, 0.40, 0.50, 0.65, 0.85, 1.20, 1.50, 2.00, 2.70, 3.50, 4.50, 6.50, 8.50, 10, 12, 15, 18.00, 22.00, 27.00, 35.00, 45.00, 60.00, 70.00, 90.00, 120.00, 150.00, 200.00, 250.00, 300.00, 400.00)#, 500.00, 650.00, 800.00, 1000.00, 1200.00, 1500.00, 2000.00, 3000.00)

# this is to cache the results
#' @export
findU1modesForAllQ2 <- function(key = NULL, A0 = NULL, h = 0.01, alpha = 0, doPlot = FALSE) {
  if(is.null(A0))
    A0 <- get('A0', envir = ihqcdEnv)

  if(is.null(key))
    key <- paste('u1modes', A0, h, alpha, sep = '-')

  # cat('\n Finding modes for A0 =', A0, ', h =', h, ', Q2s', Q2s, '\n')
  # set the constants in the IHQCD package
  # here a large zmax is used since we have been experiencing some boundary effect
  ihqcd <- iHQCD(A0 = A0, h = 0.001, zmax = 20)
  #populate the ihqcd environment
  ihqcd$solve()

  # this is the set of different Qs we have right now
  modes.names <- c("Q2", "x", "y")
  modes <- vector("list", length(modes.names))
  names(modes) <- modes.names

  if(doPlot)
    plot.new()

	cl <- rainbow(length(Q2s))
  i <- 1
  for(Q2 in Q2s) {
    m <- force(findU1Modes(Q2, alpha = alpha))

    # DEBUG
    # do not clip text outside the graph
    if(doPlot) {

      zmax <- get('zmax', envir = ihqcdEnv)
      par(xpd=NA)
      if(i == 1)
        plot(m$x, m$y, type = 'l', xlim = c(0, 1.2 * zmax), xlab = "z", ylab = "f(z)", main = "Non-normalizable modes for all Q")
      else {
        lines(m$x, m$y, type = 'l', xlim = c(0, 1.2 * zmax), col= cl[i])
        halfIndex <- (length(Q2s) - i + 1)^2 * as.integer(length(m$x) / 1) / length(Q2s)^2
        posy <- m$y[halfIndex]
        posx <- m$x[halfIndex]
        if(i > 12) {
          posy = i * 0.07 - 1
          posx = posx - 1
        }
        text(posx, posy, labels = bquote(Q^2 ==.(Q2)), col= cl[i])
      }
    }
    # END DEBUG

    modes$Q2[[i]] = Q2
    modes$x[[i]] = m$x
    modes$y[[i]] = m$y
    i <- i + 1
  }
  rredis::redisSet(key, modes);

  # set everything back
  ihqcd <- iHQCD(A0 = A0, h = h)
  #populate the ihqcd environment
  ihqcd$solve()
}

#' @export
computeU1DeltaFun <- function(Q2, A0 = NULL, h = 0.01) {
  m       <- getMode(Q2, A0, h)
  fspline <- splinefun(m$x, m$y)
  df      <- fspline(m$x, deriv = 1)

  A  <- get('A', envir = ihqcdEnv)
  As <- get('As', envir = ihqcdEnv)
  z  <- get('z', envir = ihqcdEnv)
  la <- get('lambda', envir = ihqcdEnv)
  Afun  <- splinefun(z, A)
  Asfun <- splinefun(z, As)
  lafun <- splinefun(z, la)
  # Alfonso's note eq 2.48
  t1  <- m$y^2 + (1/Q2) * df^2
  p13 <- t1 * exp(Asfun(m$x)) / lafun(m$x)
  p13  <- p13 * Q2 * exp(-2 * Afun(m$x))
  list(x = m$x, y = p13, t1 = t1, df = df)
}
