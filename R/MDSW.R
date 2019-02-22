#' @import MDSW

#DEPRECATED
modswEnv <- new.env()

#' Constructor
#' @export
MDSW <- function(muG, mq, sigma) {

  # Solves the Modified Sof-Wall model of reference arXiv 1303.6929
  # Finds the warp factor in the string frame As
  SolveMDSW  <- function(muG = 0.1157, mq = 0.0, sigma = 0.02349) {

    flog.debug(paste('[MDSW] Solving HVQCD for muG', muG, ', mq', mq, ', sigma', simga^3))
    # solve the problem and then store the relevant variables
    sol          <- solveMDSW(muG, mq, sigma)
    z            <- sol$z
    As           <- sol$As
    # compute some extras, derivatives are useful for later calculations
    Asder1        <- sol$dAs
    Asder2        <- sol$d2As
    AsSpline      <- splinefun(z, As)
    Phi           <- sol$Phi
    Phider1       <- sol$dPhi
    Phider2       <- sol$d2Phi
    Phifun        <- splinefun(z, Phi)
    Chi           <- sol$chi
    Chider1       <- sol$dchi
    Chifun        <- splinefun(z, Chi)
    # Warp factor in the Einstein frame
    A             <- As - 2 * Phi / 3.0
    Ader1         <- sol$dAs - 2 * sol$dPhi / 3.0
    Ader2         <- sol$d2As - 2 * sol$d2Phi / 3.0
    Aspline       <- splinefun(z, A)
    # compute the vector meson potential
    urho         <- 0.5 * (Asder2 - Phider2) + 0.25 * (Asder1 - Phider1) ^ 2
    # compute the scalar glueball potential
    u0           <- 0.5 * (3 * Asder2 - 2 * Phider2) + 0.25 * (3 * Asder1 - 2 * Phider1) ^ 2
    # compute the tensor glueball potential
    u2           <- 0.5 * (3 * Asder2 - 2 * Phider2) + 0.25 * (3 * Asder1 - 2 * Phider1) ^ 2 + 4 * exp(2 * As)
    # some combinations that appear in the potential of the graviton trajectory
    aF           <- Phider2
    bF           <- Asder2 - Asder1^2
    cF           <- Phider1^2
    l1_2         <- sqrt(exp(Phi))
    e2As         <- exp(2 * As)
    e2A          <- exp(2 * A)

    # to avoid z=0 boundary effects the first values are removed
    # the biggest problem comes from the derivative of A
    len <- length(z)
    st  <- 6
    list(z = z[st:len],
         u0 = u0[st:len],
         u2 = u2[st:len],
         urho = urho[st:len],
         A = A[st:len],
         Ader1 = Ader1[st:len],
         Ader2 = Ader2[st:len],
         Afun  = Aspline,
         As = As[st:len],
         Asder1 = Asder1[st:len],
         Asder2 = Asder2[st:len],
         AsSpline = AsSpline,
         Chi    = Chi[st:len],
         Chider1 = Chider1[st:len],
         Chufun = Chifun,
         Phi = Phi[st:len],
         Phider1 = Phider1[st:len],
         Phider2 = Phider2[st:len],
         Phifun  = Phifun,
         aF = aF[st:len],
         bF = bF[st:len],
         cF = cF[st:len],
         l1_2 = l1_2[st:len],
         e2As = e2As[st:len],
         e2A = e2A[st:len],
         muG = muG,
         mq  = mq,
         sigma = sigma)
  }

  loadGlobally <- function(s) {
    #flog.debug('Loading ModSW results in the global environment')
    # put all in the global environment
    mapply(function(n, v) assign(n, v, envir = .GlobalEnv), names(s), s)
  }

  i <- list(solve = cache(SolveMDSW),
            loadGlobally = loadGlobally)
  class(i) <- append(class(i), 'MDSW')
  i
}

#' Implements the generic function solves for the MDSW model
#' allowing to cache the results.
#' @export
solve.MDSW <- function(modsw, muG = 0.1157, mq = 0.0, sigma = 0.02349) {
  modsw$loadGlobally(modsw$solve(muG = muG, mq = mq, sigma = sigma))
}
