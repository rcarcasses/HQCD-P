#' @export
tachyon <- function(x = 0) {
  tauVsA<- NULL
  A      <- get('A', envir = ihqcdEnv)
  Ader1  <- get('Ader1', envir = ihqcdEnv)
  lambda <- get('lambda', envir = ihqcdEnv)
  z      <- get('z', envir = ihqcdEnv)
  # create a function that given an A gives the correspondent lambda
  lambdaFun <- splinefun(A, lambda)
  qFun      <- splinefun(A, exp(A) / Ader1)
  Cte       <-  81 * 3 ^(5/6) * 115^(4/3) * 11 / (812944 * 2^(1/6))
  # this is obtained from the asymptotic behavior of tau
  T0        <- function(tIR) tIR * exp(-max(z) * Cte)

  V0 <- 12
  V1 <- 11/(27*pi^2)
  V2 <- 4619/(46656*pi^4)
  lambda01 <- 8*pi^2
  W0 <- 3/11
  W1 <- (24 + (11-2*x)*W0)/(27*W0*pi*pi)
  W2 <- (24*(857-46*x)+(4619-1714*x+92*x*x)*W0)/(46656*pi^4*W0)
  a0 <- (12-x*W0)/8
  a1 <- (115-16*x)/(216*pi*pi)

  k       <- function(l)             1/(1 + 3 * a1 * l / 4)^(4/3)
  w       <- function(l)             k(l)
  dk      <- function(l)             -a1/((1+3*a1*l/4)^(7/3))
  Vg      <- function(l)             V0 * (1 + V1 * l + V2 * (l^2) * sqrt(1 + log(1 + l/lambda01))/(1 + l/lambda01)^(2/3))
  dVg     <- function(l)             V0*(V2*l*(11*l+12*lambda01)+4*V2*l*(2*l+3*lambda01)*log(1+l/lambda01)+6*V1*lambda01*(1+l/lambda01)^(5/3)*sqrt(1+log(1+l/lambda01)))/(6*lambda01*(1+l/lambda01)^(5/3)*sqrt(1+log(1+l/lambda01)))
  Vf      <- function(l, t)          W0 * (1 + W1 * l + W2 * (l^2)) * exp(-a0 * t^2)
  dVfl    <- function(l, t)          W0 * (W1 + 2 * W2 * l) * exp(-a0 * t^2)
  dVft    <- function(l, t)          -2 * a0 * t * Vf(l, t)
  dLnhVf1 <- function(l, t)          -0.5 * a1 / (1 + 3 * a1 * l / 4)  + (W1 + 2 * W2 * l) / (1 + W1 * l + W2 * (l^2))
  dLnhVf2 <- function(l, t)          -2 * a0 * t
  dLnhVf3 <- function(l, t)          -a1 / (1 + 3 * a1 * l / 4) + (W1 + 2 * W2 * l) / (1 + W1 * l + W2 * (l^2))
  # differential equations
  # X <- dq Y <- dlambda Z <- dtau
  BG2 <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      # cache some results
      lambda  <- lambdaFun(Time)
      Y       <- lambdaFun(Time, deriv = 1)
      q       <- qFun(Time)
      dq      <- qFun(Time, deriv = 1)
      kval    <- k(lambda)
      dkval   <- dk(lambda)
      vf      <- Vf(lambda, tau)
      dvft    <- dVft(lambda, tau)
      dlnhvf1 <- dLnhVf1(lambda, tau)
      dlnhvf2 <- dLnhVf2(lambda, tau)
      dlnhvf3 <- dLnhVf3(lambda, tau)

      # now compute the iteration step logic
      dtau <- Z
      dZ   <- -4 * Z + dq * Z / q - (kval * Z^3 / q^2) * (4 + Y * dlnhvf1)
              + Z^2 * dlnhvf2 - Z * Y * dlnhvf3 + q^2 * dvft / (kval * vf)
      val <- list(c(dtau, dZ))
      val
    })
  }

  plotq <- function(...) {
    plot(ihqcd$z, qFun(ihqcd$A), xlab = 'z', ylab = 'q', ...)
  }

  # find the tachyon profile given the lambda and A
  solve <- function(tIR = 10) {
    # get the IR boundary values for the warp factor and its derivative
    AUV     <- A[1]
    AIR     <- A[length(A)]
    Ader1IR <- Ader1[length(Ader1)]
    dtIR    <- Cte * tIR / Ader1IR
    initCon <- c(tau = tIR, Z = dtIR)
    Avalues <- seq(from = AIR, to = AUV, len = 100)
    sol <<- ode(func = BG2, y = initCon, times = Avalues, parms = NULL)
    plot(sol[, 1], sol[, 2])
    tauFunVsA <- splinefun(sol[, 1], sol[, 2])

    tauVsA <<- tauFunVsA(A)
  }

  plotTau <- function(...) {
    plot(ihqcd$z, tau, ylab = 'tau', xlab = 'z',...)
  }

  t <- list(
    solve = solve,
    plot = plotTau,
    tauVsA = function() tauVsA,
    plotq = plotq
  )
  class(t) <- append(class(t), 'tachyon')
  t
}


#' @export
plot.tachyon <- function(t, ...) t$plot(...)
