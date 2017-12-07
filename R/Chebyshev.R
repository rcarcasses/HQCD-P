#' @export
chebEnv <- new.env()
#' @export
ChebyshevSpectrum <- function(pot = list(x = seq(-15, 15, len = 1000), y = seq(-15, 15, len = 1000)^2), N = 40, numWFs = 10) {
  # map the problem to the range [-1,1]
  a <- min(pot$x)
  b <- max(pot$x)
  scal <- (b - a) / 2
  scal <- scal^2
  #cat('Interval [', a, ',', b, '] mapped to [-1,1], scale factor', scal, '\n')

  VfunUnscal <- splinefun(pot$x, pot$y)
  Vfun <- function(x) {
    scal * VfunUnscal(0.5 * ((b - a) * x + b + a))
  }
  # create the grid
  x <- -cos(pi  * seq(0, N) / N)
  L <- length(x)

  key    <- paste0('Tm', N)
  keyInt <- paste0('TmInt', N)
  if(!exists(key, envir = chebEnv)) {
    cat('Computing Tm and TmInt\n')
    # checked in mathematica with the right results
    Tm <- matrix(nrow = L, ncol = L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i == 1)      # first
          Tm[i, j] <- 1
        else if(i == 2) # second
          Tm[i, j] <- x[j]
        else
          Tm[i, j] <- 2 * x[j] * Tm[i - 1, j] - Tm[i - 2, j]
      }
    }
    # now the integral, last column turns out to be always zero
    TmInt <- matrix(nrow = L, ncol = L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i == 1)      # first
          TmInt[i, j] <- x[j] + 1
        else if(i == 2) # second
          TmInt[i, j] <- 0.5 * (x[j]^2 - 1)
        else if(i != L)
          TmInt[i, j] <- Tm[i + 1, j] / (2 * (i)) - Tm[i - 1, j] / (2 * (i - 2)) + ((-1)^i) / ((i - 1)^2 - 1)
        else   # case i = L
          TmInt[i, j] <- (2 * x[j] * Tm[i, j] - Tm[i - 1, j]) / (2 * (i)) - Tm[i - 1, j] / (2 * (i - 2)) + ((-1)^i) / ((i - 1)^2 - 1)
      }
    }

    assign(keyInt, TmInt, envir = chebEnv)
    assign(key, Tm, envir = chebEnv)
  } else {
    TmInt <- get(keyInt, envir = chebEnv)
    Tm    <- get(key, envir = chebEnv)
  }

  # compute D_ik matrix, checked vs paper
  D <- matrix(nrow = L, ncol = L)
  for (i in 1:L) {        # loop through rows
    for (k in 1:L) {      # loop through columns
      v <- unlist(sapply(seq(1, L), function(j) {   # primed '' sum
                                              f <- 1
                                              if(j == 1 || j == L)
                                                f <- 0.5
                                              f * Tm[j, k] * TmInt[j, i]
                                          }))
      D[i, k] <- (2 / N) * sum(v)
      # prime '' sum
      if(k == 1 || k == L)
        D[i, k] <- 0.5 * D[i, k]
    }
  }

  u      <- 0.5 * (x + 1)
  E      <- D %*% D

  etilde <- E[L,]
  V      <- sapply(x, Vfun)
  Ehat   <- matrix(nrow = L, ncol = L)
  Dhat   <- matrix(nrow = L, ncol = L)
  for (i in 1:L)         # loop through rows
    for (j in 1:L) {       # loop through columns
      Ehat[i, j] <- E[i, j] * V[j]
      Dhat[i, j] <- D[i, j] * V[j]
    }

  S      <- D %*% Dhat
  #actually S and Ehat are the same as expected
  # cat('values\n')
  # print(S - Ehat)
  stilde <- S[L,]

  US     <- matrix(nrow = L, ncol = L)
  UE     <- matrix(nrow = L, ncol = L)
  for (i in 1:L)         # loop through rows
    for (j in 1:L) {       # loop through columns
      UE[i, j]   <- u[i] * etilde[j]
      US[i, j]   <- u[i] * stilde[j]
    }

  #cat('\nUS = \n', US)
  # compute the A and B matrices
  Id <- diag(1, nrow = L, ncol = L)
  A <- -Id - US + Ehat
  B <- E - UE
  # now we remove the first and last row/columns since in B they are 0
  B   <- B[2:(L - 1), 2:(L - 1)]
  A   <- A[2:(L - 1), 2:(L - 1)]
  eg  <- eigen(solve(B) %*% A)
  wfs <- eg$vectors[,(L - 2) : (L - 1 - numWFs)]
  # finally we need to translate back the values of the wavefunctions to the original grid
  xval <- 0.5 * ((b - a) * x + b + a)
  xval <- xval[2:(L - 1)]
  wfs <- lapply(1:numWFs, function(i) list(x = xval, y = wfs[,i]))
  # solve the eigenvalue proble
  cb <- list(energies = eg$values[(L - 2) : (L - 1 - numWFs)]/ scal,
             wfs = wfs)
  class(cb) <- append(class(cb), 'ChebyshevSpectrum')
  cb
}
