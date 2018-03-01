odeFunAdS <- function(x, y, pars) {
  dy1 <- y[2]
  dy2 <- Q2 * y[1] + y[2] / x
  list(c(dy1, dy2))
}

odeFun1 <- function(x, y, pars) {
  dy2 <- Q2 * y[1] - (Asder1fun(x) - Phider1fun(x)) * y[2]
  dy1 <- y[2]
  list(c(dy1, dy2))
}

odeFun2 <- function(x, y, pars) {
  # remark: we are using Q2=q^2 and not the usual -q^2 since we are using the GR metric convention for eta
  dy2 <- Q2 * y[1] * fact0fun(x) - (Asder1fun(x) - Phider1fun(x)
            + 2 * alpha * (2 * Asder1fun(x) * Asder2fun(x) - Asder3fun(x)) / fact1fun(x)) * y[2]
  dy1 <- y[2]
  list(c(dy1, dy2))
}

postcompute <- function(pre) {
  # use the default solution globally again
  solve(pre$ihqcd, A0 = pre$originalPars$A0, h = pre$originalPars$h, zmax = pre$originalPars$zmax)
}

precompute <- function(Q2, alpha) {
  ihqcd <- iHQCD()
  originalPars <- list(A0 = A0, h = h, zmax = zmax)
  solve(ihqcd, A0 = A0, h = 0.001, zmax = 20)
  Asfun        <- splinefun(z, As)
  Asder1fun    <- splinefun(z, Asder1)
  Asder2fun    <- splinefun(z, Asder2)
  Asder3fun    <- splinefun(z, Asder3)
  Phider1fun   <- splinefun(z, Phider1)
  fact0fun     <- splinefun(z, (1 - 2 * exp(-2 * As) * alpha * Asder2) / (1 - 2 * exp(-2 * As) * alpha * Asder1^2))
  fact1fun     <- splinefun(z, exp(2 * As) - 2 * alpha * Asder2)
  list(ihqcd = ihqcd, originalPars = originalPars , env = environment())
}

# these are the non-normalizable modes
getU1NNModeRaw <- function(Q2 = 3.5, alpha = 0) {
  flog.debug(paste('finding mode Q2', Q2, 'alpha', alpha))
  # first we need to compute ihqcd with bigger precision
  # set the constants in the IHQCD package
  # here a large zmax is used since we have been experiencing some boundary effect
  pre <- precompute(Q2, alpha)
  # this is useful for the non minimal coupling case
  # now we need to define the differential equation for the U(1) field modes
  odeFun <- odeFun1
  if (abs(alpha) > 0)
    odeFun <- odeFun2

  # pass all the variables inside the environment of this function
  environment(odeFun) <- pre$env
  # get the solution as a data frame
  m <- as.data.frame(bvpcol(yini = c(y = 1, dy = NA), x = z, func = odeFun, yend = c(0, NA), nmax = 100000, atol = 1e-9))
  mode <- getFQs(Q2, m)
  # roll back global changes
  postcompute(pre)
  mode
}

getFQs <- function(Q2, m) {
  # also compute some useful combinations
  fQ      <- splinefun(m$x, m$y)
  fQ2     <- splinefun(m$x, m$y^2)
  dfQ     <- splinefun(m$x, m$dy)
  dfQ2    <- splinefun(m$x, m$dy^2)
  factor  <- splinefun(m$x, m$y^2 + (1/Q2) * m$dy^2)
  list(fQ = fQ, fQ2 = fQ2, dfQ2 = dfQ2, factor = factor)
}

getFQNormalizableModes <- function(mass, m) {
  # also compute some useful combinations
  fmsplin <- splinefun(m$z, m$f)
  fm      <- Vectorize(function(x) if (x <= m$z[length(m$z)]) fmsplin(x) else 0)
  dfm     <- Vectorize(function(x) if (x <= m$z[length(m$z)]) fmsplin(x, deriv = 1) else 0)
  list(fm = fm, dfm = dfm)
}

# these are the normalizable modes
getU1NormalizableModeRaw <- function(mass = 1.0, alpha = 0) {
  flog.debug(paste('finding mode mass', mass, 'alpha', alpha))
  # first we need to compute ihqcd with bigger precision
  # set the constants in the IHQCD package
  # here a large zmax is used since we have been experiencing some boundary effect
  pre <- precompute(-mass^2, alpha)
  # this is useful for the non minimal coupling case
  # now we need to define the differential equation for the U(1) field modes
  odeFun <- odeFun1

  # TODO - CHECK THAT IF ODEFUN2 NEEDS TO BE MODIFIED
  if (abs(alpha) > 0)
    odeFun <- odeFun2

  # pass all the variables inside the environment of this function
  environment(odeFun) <- pre$env
  # Set up parameters and boundary conditions
  pars <- c(mass = mass, alpha = alpha)
  # In the UV warped AdS is the same as AdS
  # So we use AdS solution as boundary condition
  y0AdS    <- mass * z[1] * besselJ(mass * z[1], 1)
  dy0AdS   <- mass * mass * z[1] * besselJ(mass * z[1], 0)
  yini <- c(y1 = y0AdS, y2 = dy0AdS)
  # get the solution as a data frame
  m <- as.data.frame(ode(y = yini, times = z, func = odeFun, parms = pars))
  # Now we need to normalize the normalizable mode
  # We do that by imposing df = 0 at some zf
  f <- splinefun(m$time, m$y1)
  # We now the value of zf that satisfied the BC condition
  df <- splinefun(m$time, m$y2)
  zf <- min(uniroot.all(df,
                lower = m$time[1], upper = 2.4048 * 1.5/ mass))
  z    <- m$time[m$time < zf]
  y1   <- m$y1[m$time < zf]
  y2   <- m$y2[m$time < zf]
  ff <- function(z) exp(Asder1fun(z) - Phider1fun(z))
  environment(ff) <- pre$env
  f <- splinefun(z, ff(z) * y1 * y1)
  C2 <- integrate(f = f, lower = z[1], upper = z[length(z)])$value
  f <- y1 / sqrt(C2)
  df <- y2 / sqrt(C2)
  m <- data.frame(z, f, df)
  mode <- getFQNormalizableModes(mass, m)
  # roll back global changes
  postcompute(pre)
  # This solution does not have the correct
  mode
}

#' @export
getU1NNMode <- cache(getU1NNModeRaw)

#' @export
getU1NormalizableMode <- cache(getU1NormalizableModeRaw)
