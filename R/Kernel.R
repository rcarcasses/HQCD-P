# kernelName is important when doing multikernel computations, should be unique
# numReg indicates the number of reggions to be considered while finding the kernel
#' @export
kernelBCDQ <- function(potPars = 'b', numReg = 3, kernelName = '') {
  p <- potential(potPars)

  pot <- p$u()
  z   <- p$z
  lastKernel <- NULL
  time.taken <- NULL

  # set the potential to an external one
  setNewPotential <- function(f, pars) {
    p$setNewPotential(f, pars)
    # update the local potential function
    pot <<- p$u()
  }

  # the arguments here should be consistent with the potPars passed
  findKernelFun <- function(...) {
    fArgs <- do.call(p$extractPotentialParameters, list(...))    # sometimes the coefficients are passed to this, just ignore them
    start.time <- Sys.time()
    t <- function(J, n) {
      args <- as.list(c(J = J, fArgs))
      u2j  <- do.call(pot, args)
      dE   <- abs(min(u2j)) / 400;
      # make it at least one for values that give a minimum near 0
      dE   <- max(dE, 0.01)
      # cat("j =", J, " d0 =", dE, "\n")
      setPotential(z, u2j)
      computeSpectrum(n, dE)
      energies <- getEnergies()
      wf <- getWavefunctions()[[n]]
      list(t = energies[[n]], wf = wf, u2j = u2j, energies = energies)
    }

    tvec <- function(J, n) sapply(J, function(j) t(j, n)$t)

    # define the get intecept function
    Js <- seq(0.2, 2, len = 20)
    getIntercept <- function(n) {
      tspline <- function(j) tvec(j, n)
      roots <- uniroot.all(tspline, c(0.2, 2.2), n = 6)
      if(length(roots) > 0)
        js <- max(roots)
      else
        js <- 0

      tresult <- t(js, n)
      list(js = js, wf = tresult$wf, u2j = tresult$u2j)
    }

    if(!is.null(cl))
      s <- lapply(seq(1, numReg), getIntercept)
    else
      s <- parLapply(cl, seq(1, numReg), getIntercept)

    # normalize the wave function
    lapply(1:numReg, function(i){
      # cat(paste('j', i, sep = ''),'= ', s$js, '     ')
      wffun <- splinefun(s[[i]]$wf$x, s[[i]]$wf$y^2)
      c <- integrate(wffun, s[[i]]$wf$x[1], s[[i]]$wf$x[length(s[[i]]$wf$x)])$value
      s[[i]]$wf$y <<- (s[[i]]$wf$y[10] / abs(s[[i]]$wf$y[10])) * s[[i]]$wf$y / sqrt(c)
    })

    end.time <- Sys.time()
    time.taken <<- end.time - start.time

    # we also want to compute the spectral data to compare with the glueball spectrum (or mesons)
    # so far we need only to compare with the first two levels
    gs  <- lapply(2 * (1:3), function(J) t(J, 2))
    tr1 <- lapply(gs, function(v) v$energies[[1]])
    tr2 <- lapply(gs, function(v) v$energies[[2]])

    lastKernel <<- as.list(c(s, list(tr1 = tr1, tr2 = tr2)))
    lastKernel
  }

  findKernel <- rediscache::cache(findKernelFun, kernelName)
  extractPotentialParameters <- function(...) do.call(p$extractPotentialParameters, list(...))

  k <- list(findKernel = findKernel,
            potential = p,
            setNewPotential = setNewPotential,
            extractPotentialParameters = extractPotentialParameters,
            getLast = function() lastKernel,
            getTimeTaken = function() time.taken)
  class(k) <- append(class(k), 'kernelBCDQ')
  k
}
