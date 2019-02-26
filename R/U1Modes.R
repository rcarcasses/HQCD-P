odeFunAdS <- function(x, y, pars) {
  dy1 <- y[2]
  dy2 <- Q2 * y[1] + y[2] / x
  list(c(dy1, dy2))
}

odeFun1 <- function(x, y, pars) {
  dy2 <- Q2 * y[1] - (AsSpline(x, deriv = 1) - Phifun(x, deriv = 1)) * y[2]
  dy1 <- y[2]
  list(c(dy1, dy2))
}

odeFun2 <- function(x, y, pars) {
  # remark: we are using Q2=q^2 and not the usual -q^2 since we are using the GR metric convention for eta
  dy2 <- Q2 * y[1] * fact0fun(x) - (AsSpline(x, deriv = 1) - Phifun(x, deriv = 1)
            + 2 * alpha * (2 * AsSpline(x, deriv = 1) * AsSpline(x, deriv = 2) - AsSpline(x, deriv = 3)) / fact1fun(x)) * y[2]
  dy1 <- y[2]
  list(c(dy1, dy2))
}

postcompute <- function(pre) {
  if (model == "iHQCD")
  {
    # use the default solution globally again
    solve(pre$ihqcd, A0 = pre$originalPars$A0, h = pre$originalPars$h, zmax = pre$originalPars$zmax)
  }
}

precompute <- function(Q2, alpha) {
  if (model == "iHQCD")
  {
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
  else list(env = environment())
}

# these are the non-normalizable modes
getU1NNModeRaw <- function(Q2 = 3.5, alpha = 0) {
  flog.debug(paste('finding mode Q2', Q2, 'alpha', alpha))
  # In the case of iHQCD we need to compute iHQCD with bigger precision
  # set the constants in the IHQCD package
  # here a large zmax is used since we have been experiencing some boundary effect
  # this is useful for the non minimal coupling case
  pre <- precompute(Q2, alpha)
  # we define the differential equation for the U(1) field modes
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

getFQNormalizableModes <- function(m) {
  # also compute some useful combinations
  fm  <- splinefun(m$z, m$f)
  dfm <- splinefun(m$z, fm(m$z, deriv = 1) )
  return(list(fm = fm, dfm = dfm))
}

# these are the normalizable modes
getU1NormalizableModeRaw <- function()
{
  if (model == "MDSW")
  {
    s <- computeSpectrum(z, urho, nEigen = 1)
    # Get the xx coordinates
    xx <- s$wfs[[1]]$x
    # Get the yy coordinates
    yy <- s$wfs[[1]]$y
    # Computation of the interpolation function
    psifun <- splinefun(xx, yy)
    integral <- integrate(function(x) psifun(x) * psifun(x), z[1], z[length(z)])$value
    yy <- yy / sqrt(integral)
    m <- list(z = xx, f = yy)
    mode <- getFQNormalizableModes(m)
    return(mode)
  }
  else print("Normalizable modes not defined for this model.")
}


#' @export
getU1NNMode <- cache(getU1NNModeRaw)

#' @export
getU1NormalizableMode <- cache(getU1NormalizableModeRaw)
