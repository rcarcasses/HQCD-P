#' @import HVQCD

#DEPRECATED
hvqcdEnv <- new.env()

#' Constructor
#' @export
HVQCD <- function(x = 1.0, t0 = 1.0) {
  # Solves the Kiritsis and Jarvinen Model arXiv:1112.1261
  # Finds z, dz, lambda, tau, dlambda and dtau as a function of A. mq is also returned
  SolveHVQCD  <- function(x = 1.0, t0 = 1.0) {

    flog.debug(paste('[HVQCD] Solving HVQCD for x', x, ', t0', t0))
    # solve the problem and then store the relevant variables
    sol <- solveHVQCD(x, t0)
    # Store the values of z, A, lambda and tau
    z <- sol$z
    A <- sol$A
    dz <- sol$dz
    lambda <- sol$lambda
    tau <- sol$tau
    mq <- sol$mq
    # compute some extras
    Aspline      <- splinefun(z, A)
    lambdaspline <- splinefun(z, lambda)
    tauspline    <- splinefun(z, tau)
    Phi          <- log(sol$lambda)
    Phispline    <- splinefun(z, Phi)
    tauspline    <- splinefun(z, tau)
    # derivatives are useful for different applications, so let's export them as well
    Ader1        <- Aspline(z, deriv = 1)
    Ader2        <- Aspline(z, deriv = 2)
    Ader3        <- Aspline(z, deriv = 3)
    lambdader1   <- lambdaspline(z, deriv = 1)
    lambdader2   <- lambdaspline(z, deriv = 2)
    Phider1      <- Phispline(z, deriv = 1)
    Phider2      <- Phispline(z, deriv = 2)
    tauder1      <- tauspline(z, deriv = 1)
    tauder2      <- tauspline(z, deriv = 2)
    # also sometimes the warp factor in the string frame is needed
    As           <- A + (2/3) * Phi
    AsSpline     <- splinefun(z, As)
    Asder1       <- AsSpline(z, deriv = 1)
    Asder2       <- AsSpline(z, deriv = 2)
    Asder3       <- AsSpline(z, deriv = 3)
    st  <- 10
    len <- length(z)-st
    list(z = z[st:len],
         A = A[st:len],
         Ader1 = Ader1[st:len],
         Ader2 = Ader2[st:len],
         Ader3 = Ader3[st:len],
         Afun  = Aspline,
         lambda = lambda[st:len],
         lambdader1 = lambdader1[st:len],
         lambdader2 = lambdader2[st:len],
         Phi = Phi[st:len],
         Phider1 = Phider1[st:len],
         Phider2 = Phider2[st:len],
         Phifun  = Phispline,
         tau     = tau[st:len],
         tauder1 = tauder1[st:len],
         tauder2 = tauder2[st:len],
         taufun  = tauspline,
         As = As[st:len],
         AsSpline = AsSpline,
         Asder1 = Asder1[st:len],
         Asder2 = Asder2[st:len],
         Asder3 = Asder3[st:len],
         mq = mq)
  }

  loadGlobally <- function(s) {
    #flog.debug('Loading HVQCD results in the global environment')
    # put all in the global environment
    mapply(function(n, v) assign(n, v, envir = .GlobalEnv), names(s), s)
  }

  i <- list(solve = cache(SolveHVQCD),
            loadGlobally = loadGlobally)
  class(i) <- append(class(i), 'HVQCD')
  i
}

#' Implements the generic function solves for the HVQCD model
#' allowing to cache the results.
#' @export
solve.HVQCD <- function(hvqcd, x = 1.0, t0 = 1.0) {
  hvqcd$loadGlobally(hvqcd$solve(x = x, t0 = t0))
}
