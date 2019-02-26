#' @import HVQCD

#DEPRECATED
hvqcdEnv <- new.env()

#' Constructor
#' @export
HVQCD <- function(x = 1.0, t0 = 1.0, W0 = 12.0/11,  V0 = 12, lambda0 = 8 * pi^2) {
  # Solves the Kiritsis and Jarvinen Model arXiv:1112.1261
  # Finds z, dz, lambda, tau, dlambda and dtau as a function of A. mq is also returned
  SolveHVQCD  <- function(x = 1.0, t0 = 1.0, W0 = 12.0/11,  V0 = 12, lambda0 = 8 * pi^2) {

    flog.debug(paste('[HVQCD] Solving HVQCD for x', x, ', t0', t0, ', W0', W0, ', V0', V0, ', lambda0', lambda0))
    # solve the problem and then store the relevant variables
    sol <- solveHVQCD(x, t0, W0, V0, lambda0)
    # Store the values of z, A, lambda and tau
    z <- sol$z
    A <- sol$A
    lambda <- sol$lambda
    tau <- sol$tau
    mq <- sol$mq
    # compute some extras
    Aspline      <- splinefun(z, A)
    lambdaspline <- splinefun(z, lambda)
    tauspline    <- splinefun(z, tau)
    Phi          <- log(lambda)
    Phispline    <- splinefun(z, Phi)
    # derivatives are useful for different applications, so let's export them as well
    Ader1        <- Aspline(z, deriv = 1)
    Ader2        <- Aspline(z, deriv = 2)
    Ader3        <- Aspline(z, deriv = 3)
    lambdader1   <- lambdaspline(z, deriv = 1)
    lambdader2   <- lambdaspline(z, deriv = 2)
    Phider1      <- Phispline(z, deriv = 1)
    Phider2      <- Phispline(z, deriv = 2)
    Phider3      <- Phispline(z, deriv = 3)
    tauder1      <- tauspline(z, deriv = 1)
    tauder2      <- tauspline(z, deriv = 2)
    # also sometimes the warp factor in the string frame is needed
    As           <- A + (2/3) * Phi
    AsSpline     <- splinefun(z, As)
    Asder1       <- AsSpline(z, deriv = 1)
    Asder2       <- AsSpline(z, deriv = 2)
    Asder3       <- AsSpline(z, deriv = 3)
    # compute the scalar glueball potential
    u0           <- (9/4) * Ader1^2 - (3/2) * Ader2 + 2 * (Ader2/Ader1)^2 +
                    3 * Ader1 * Phider2 / Phider1 - 2 * Ader2 * Phider2 / (Ader1 * Phider1)
                    - Ader3 / Ader1 + Phider3 / Phider1
    # compute the tensor glueball potential
    u2           <- (3/2) * Ader2 + (9/4) * Ader1^2
    # Quantities that are useful to compute potentials of the vector mesons
    V0 <- 12
    V1 <- 11/(27*pi^2)
    V2 <- 4619/(46656*pi^4)
    lambda0 <- 8*pi^2
    W0 <- 12/11
    W1 <- (24 + (11-2*x)*W0)/(27*W0*pi^2)
    W2 <- (24*(857-46*x)+(4619-1714*x+92*x^2)*W0)/(46656*pi^4*W0)
    a0 <- (12-x*W0)/8
    a1 <- (115-16*x)/(216*pi^2)
    # Values of the potentials useful to later compute the Vector Meson Spectra
    k  <- 1/(1 + 3 * a1 * lambda / 4)^(4/3)
    Vf <- W0 * (1 + W1 * lambda + W2 * (lambda^2)) * exp(-a0 * tau^2)
    # Function needed to compute the Vector Meson Spectra
    G            <- sqrt(1 + exp(-2 * A) * k * tauder1^2)
    Gspline      <- splinefun(z, G)
    dG           <- Gspline(z, deriv = 1)
    # some combinations that appear in the potential
    aF           <- Phider2
    bF           <- Asder2 - Asder1^2
    cF           <- Phider1^2
    l1_2         <- sqrt(lambda)
    e2As         <- exp(2 * As)
    e2A          <- exp(2 * A)
    ChiV <- k * sqrt(Vf * exp(A))
    st  <- 12
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
         Phider3 = Phider3[st:len],
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
         u0     = u0[st:len],
         u2     = u2[st:len],
         aF     = aF[st:len],
         bF     = bF[st:len],
         cF     = cF[st:len],
         l1_2   = l1_2[st:len],
         e2As   = e2As[st:len],
         e2A    = e2A[st:len],
         e2As   = e2As[st:len],
         V0     = V0,
         V1     = V1,
         V2     = V2,
         lambda0 = lambda0,
         W0     = W0,
         W1     = W1,
         W2     = W2,
         a0     = a0,
         a1     = a1,
         k      = k[st:len],
         Vf     = Vf[st:len],
         G      = G[st:len],
         dG     = dG[st:len],
         Gfun   = Gspline,
         ChiV   = ChiV[st:len],
         mq = mq,
         model = "HVQCD")
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
solve.HVQCD <- function(hvqcd, x = 1.0, t0 = 1.0, V0 = 12, W0 = 12.0/11, lambda0 = 8 * pi^2) {
  hvqcd$loadGlobally(hvqcd$solve(x = x, t0 = t0, V0 = V0, W0 = W0, lambda0 = lambda0))
}
