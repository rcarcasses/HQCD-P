#' A kernel unit represents a single type of pomeron kernel
#' @param  potential A function that computes the potential and has as arguments J plus the potential parameters
#' @param  numReg Indicates the number of reggions to be considered while finding the kernel
#' @param  kernelName Is important when doing multikernel computations, should be unique
#' @param comment Can be used to add some additional information about what this kernel represents
#' @param optimPars Indicate which parameters can be optimized for this kernel
#' @export
kernelUnit <- function(potential, numReg = 3, kernelName = '', comment = '', optimPars = c()) {
  lastKernel <- NULL
  time.taken <- NULL

  # the arguments here should be consistent with the potPars passed
  # the .t name is to prevent partial matching with extra arguments passed
  # (please prevent users to use variable names starting with .!)
  findKernelFun <- function(.t = 0, ...) {
    flog.debug(paste('finding kernel for t =', .t, ' pars:', do.call(paste, as.list(format(list(...), digits = 3)))))
    # extra parameters may come, just ignore them
    fArgs <- getPotentialArgs(potential, ...)
    start.time <- Sys.time()
    t <- function(J, n) {
      u2j  <- do.call(UJgTest, as.list(c(J = J, fArgs)))
      s <- computeSpectrum(z, u2j, n)
      list(t = s$energies[[n]], wf = s$wf[[n]], u2j = u2j, energies = s$energies)
    }

    tvec <- function(J, n) sapply(J, function(j) t(j, n)$t)

    # define the get intecept function
    Js <- seq(0.2, 2, len = 20)
    getIntercept <- function(n) {
      tspline <- function(j) tvec(j, n) - .t
      roots <- uniroot.all(tspline, c(0.2, 2.2), n = 6)
      if(length(roots) > 0)
        js <- max(roots)
      else
        js <- 0

      tresult <- t(js, n)
      list(js = js, wf = tresult$wf, u2j = tresult$u2j, name = paste0(kernelName, '.', n))
    }

    # use mcapply to parallelize the computation
    # The current implementation of schrodinger package
    # does not work with mclapply in mac
    lapplyType <- lapply #if(Sys.info()['sysname'] == 'Linux') mclapply else lapply
    s <- lapplyType(1:numReg, getIntercept)

    # normalize the wave function
    lapply(1:numReg, function(i){
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

  findKernel <- cache(findKernelFun, kernelName)

  k <- list(findKernel = findKernel,
            potential = potential,
            comment = comment,
            optimPars = optimPars,
            numReg = numReg,
            getLast = function() lastKernel,
            getTimeTaken = function() time.taken)
  class(k) <- append(class(k), 'kernelUnit')
  k
}

#' @export
getPotentialArgs <- function(pot, ...) {
  fArgs <- list(...)
  fDef  <- formalArgs(pot)
  fDef  <- fDef[-match('J', fDef)]  # J is not a parameter
  fArgs <- fArgs[fDef]
  # remove the null values
  fArgs <- Filter(Negate(is.null), fArgs)
  fArgs
}
