#' @import deSolve

#' @export
ihqcdEnv <- new.env()

#' Constructor
#' @export
iHQCD <- function(A0 = 5#4.623038 #2.875115
                  , zmax = 10, h = 0.002, lambda0 = 0.0337462) {
    z     <- seq(exp(-A0), zmax, by = h)
    #z     <- seq(h, zmax, by = h)
    zmin  <- min(z)
    b0    <- 4.2
    b1    <- 51 * b0^2 /121
    alpha <- 2
    aa    <- (3/8) * (alpha - 1) / alpha
    LAdS  <- 1

    # Solves the Gursoy, Kiritsis and Nitti model.
    # Finds A and lambda as a function of z.
    solveIHQCD  <- function() {
        flog.info(paste('[IHQCD] Solving IHQCD for A0', A0,' in z = [', round(zmin, digits = 4), ',', zmax,';', h, '] grid size', length(z),'\n'))
        pars    <- c(b0 = b0, b1 = b1, aa = aa, LAdS = LAdS)
        initCon <- c(A = A0, lambda = lambda0)

        # differential equations
        BG1 <- function(Time, State, Pars) {
          with(as.list(c(State, Pars)), {
             X <- -b0 * lambda / (3 + 2 * b0 * lambda) - (2 * b0^2 + 3 * b1) * lambda^2 / (9 * (1 + lambda^2) * (1 + (2 * b0^2 + 3 * b1) * log(1 + lambda^2) / (18 * aa)))
             W <- (9/4) * (1 + (2/3) * b0 * lambda)^(2/3) * (1 + (2 * b0^2 + 3 * b1) * log(1 + lambda^2) / (18 * aa))^(4 * aa / 3)
             dA <- -(4/9) * W * exp(A) / LAdS
             dlambda <- -(4/3) * X * W * lambda * exp(A) / LAdS

             return(list(c(dA, dlambda)))
           })
        }

        # solve the problem and then store the relevant variables
        sol          <- ode(func = BG1, y = initCon, parms = pars, times = z)
        A            <- sol[,2]
        lambda       <- sol[,3]
        # write the result to the environment variable
        setBackground(z, A, lambda)
    }

    setBackground <- function(z, A, lambda, tau = NULL) {
        Aspline      <- splinefun(z, A)
        lambdaspline <- splinefun(z, lambda)
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
        # also sometimes the warp factor in the string frame is needed
        As           <- A + (2/3) * Phi
        AsSpline     <- splinefun(z, As)
        Asder1       <- AsSpline(z, deriv = 1)
        Asder2       <- AsSpline(z, deriv = 2)
        Asder3       <- AsSpline(z, deriv = 3)
        dress        <- exp(2 * A)
        # compute the tensor glueball potential
        u2           <- (3/2) * Ader2 + (9/4) * Ader1^2
        LambdaQCD    <- exp(A0 - 1 / (b0 * lambda0)) / (b0 * lambda0)^(b1 / b0^2)
        # some combinations that appear in the potential
        aF           <- Phider2
        bF           <- Asder2 - Asder1^2
        cF           <- Phider1^2
        l1_2         <- sqrt(lambda)
        e2As         <- exp(2 * As)
        e2A          <- exp(2 * A)

        # to avoid z=0 boundary effects the first values are removed
        # the biggest problem comes from the derivative of A
        len <- length(z)
        st  <- 6
        s   <- list(z = z[st:len],
                  u2 = u2[st:len],
                  A = A[st:len],
                  Ader1 = Ader1[st:len],
                  Ader2 = Ader2[st:len],
                  Ader3 = Ader2[st:len],
                  Asder1 = Asder1[st:len],
                  Asder2 = Asder2[st:len],
                  Asder3 = Asder3[st:len],
                  lambda = lambda[st:len],
                  lambdader1 = lambdader1[st:len],
                  lambdader2 = lambdader2[st:len],
                  Phi = Phi,
                  Phider1 = Phider1[st:len],
                  Phider2 = Phider2[st:len],
                  As = As[st:len],
                  dress = dress[st:len],
                  aF = aF[st:len],
                  bF = bF[st:len],
                  cF = cF[st:len],
                  l1_2 = l1_2[st:len],
                  e2As = e2As[st:len],
                  e2A = e2A[st:len],
                  h = h,
                  A0 = A0,
                  LambdaQCD = LambdaQCD)
        # set all these inside the environmet
        mapply(function(n, v) assign(n, v, envir = ihqcdEnv), names(s), s)
        s
    }

    # Returns the scalar potential
    U2 <- function(Aspline) {
        # compute the derivatives
        Ader1 <- Aspline(z, deriv = 1)
        Ader2 <- Aspline(z, deriv = 2)
        # and finally the tensor glueball schrodinger problem potential
        return((3/2) * Ader2 + (9/4) * Ader1^2)
    }

    i <- list(solve = solveIHQCD,
              setBackground = setBackground,
              get = function(attr) get(attr, envir = ihqcdEnv),
              U2 = U2)
    class(i) <- append(class(i), 'IHQCD')
    i
}

# find the best value for A0 such that the mass of the spin 2 glueball is 1.475 GeV
#' @export
findBestA0 <- function(startA0 = 2.8) {
  A0old <- 2.875115
  zMin  <- exp(-A0old)
  df <- function(pars) {
    # solve IHQCD
    s  <- do.call(iHQCD, as.list(pars))$solve()
    #u2 <- s$u2
    z  <- s$z
    Af <- splinefun(z, s$A)
    lambdaf <- splinefun(z, s$lambda)

    val <- (Af(zMin) - A0old)^2 + (lambdaf(zMin) - 0.05)^2
    cat('A0', pars[1], 'lambda', pars[2], 'val', val, '\n')
    val
    #setPotential(z, u2)
    # compute the mass of the lowest spin 2 exictation
    #computeSpectrum(2)
    #m2 <- getEnergies()[[2]]
    # return the difference
    #m2 - 2.055^2    # 2.055 is not a real glueball mass, but the value of m2 we get if we fix m0=1.475
  }  # find the best A0
  op <- optim(list(A0 = 5, lambda0 = 0.03), df, method = 'L-BFGS-B',
              lower = c(3, 0.02),
              upper = c(6, 0.05)
              )
  #roots <- uniroot.all(Vectorize(diff), c(2.2, 3), n = 6)
  #cat('Best A0', max(roots), '\n')
  op
}

#' @export
getAvailableSymbolsIHQCD <- function() names(as.list(ihqcdEnv))

#' @export
findBestA0byM0 <- function() {

  diff <- function(A0) {
    # solve IHQCD
    s  <- do.call(iHQCD, list(A0 = A0))$solve()
    u2 <- s$u2
    z  <- s$z
    setPotential(z, u2)
    # compute the mass of the lowest spin 2 exictation
    computeSpectrum(1)
    m2 <- getEnergies()[[1]]
    # return the difference
    m2 - 2.055^2    # 2.055 is not a real glueball mass, but the value of m2 we get if we fix m0=1.475
  }  # find the best A0

  roots <- uniroot.all(Vectorize(diff), c(2.2, 6), n = 6)
  cat('Best A0', max(roots), '\n')
  max(roots)
}
