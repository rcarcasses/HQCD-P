#' A kernel unit represents a single type of pomeron kernel
#' @param  potential A function that computes the potential and has as arguments J plus the potential parameters
#' @param  numReg Indicates the number of reggions to be considered while finding the kernel
#' @param  kernelName Is important when doing multikernel computations, should be unique
#' @param comment Can be used to add some additional information about what this kernel represents
#' @param optimPars Indicate which parameters can be optimized for this kernel
#' @export
kernelUnit <- function(potential, numReg = 3, kernelName = '', comment = '', optimPars = c()) {
  # the .t name is to prevent partial matching with extra arguments passed
  # (please prevent users to use variable names starting with .!)
  findReggeonDataFun <- function(.t = 0, n = 1, ...) {
    flog.trace(paste('finding reggeon ', n, ' for t =', .t, ' pars:', do.call(paste, as.list(format(list(...), digits = 3)))))
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
    tfun <- function(j) tvec(j, n) - .t
    roots <- uniroot.all(tfun, c(0.2, 2.2), n = 6)
    if(length(roots) > 0)
      js <- max(roots)
    else
      js <- 0

    s <- t(js, n)
    # create an interval around the js found
    intJs <- seq(js - 0.1, js + 0.1, len = 9)
    # now find numerically J(t)
    Jfun <- splinefun(t(intJs, n), intJs)
    # now find the derivative at the value of js found
    dJdt <- Jfun(.t, deriv = 1)
    # normalize the wave function
    wffun  <- splinefun(s$wf$x, s$wf$y^2)
    c      <- integrate(wffun, s$wf$x[1], s$wf$x[length(s$wf$x)])$value
    s$wf$y <- (s$wf$y[10] / abs(s$wf$y[10])) * s$wf$y / sqrt(c)
    # return a list with the data
    list(js = js, dJdt = dJdt, wf = s$wf, u2j = s$u2j, name = paste0(kernelName, '.', n))
  }

  findKernelFun <- function(.t = 0, n = 1, ...) {
    # use parallelization if possible
    lapplyType <- if(Sys.info()['sysname'] == 'Linux') mclapply else lapply
    s <- lapplyType(1:numReg, findReggeonDataFun, .t, n, ...)
    s
  }

  findKernel <- cache(findKernelFun, kernelName)

  k <- list(findKernel = findKernel,
            potential = potential,
            comment = comment,
            optimPars = optimPars,
            numReg = numReg)
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
