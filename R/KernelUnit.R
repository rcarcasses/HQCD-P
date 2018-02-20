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
  findReggeonDataFun <- function(n = 1, .t = 0, ...) {
    flog.trace(paste('finding reggeon ', n, ' for t =', .t, ' pars:', do.call(paste, as.list(format(list(...), digits = 3)))))
    # extra parameters may come, just ignore them
    fArgs <- getPotentialArgs(potential, ...)
    start.time <- Sys.time()
    t <- function(J, fullAns = FALSE) {
      u2j  <- do.call(UJgTest, as.list(c(J = J, fArgs)))
      s <- computeSpectrum(z, u2j, n)
      if(fullAns)
        list(t = s$energies[[n]], wf = s$wf[[n]], u2j = u2j, energies = s$energies)
      else
        s$energies[[n]]
    }

    # define the get intecept function
    tvec <- Vectorize(function(J) t(J) - .t)
    roots <- uniroot.all(tvec, c(0.2, 2.2), n = 6)
    if(length(roots) > 0)
      js <- max(roots)
    else
      js <- 0
    # now get the full answer
    s <- t(js, TRUE)
    # create an interval around the js found
    intJs <- seq(js - 0.1, js + 0.1, len = 9)
    # now find numerically J(t)
    Jfun <- splinefun(tvec(intJs), intJs)
    # now find the derivative at the value of js found
    dJdt <- Jfun(.t, deriv = 1)
    # normalize the wave function
    wffun  <- splinefun(s$wf$x, s$wf$y^2)
    c      <- integrate(wffun, s$wf$x[1], s$wf$x[length(s$wf$x)])$value
    s$wf$y <- (s$wf$y[10] / abs(s$wf$y[10])) * s$wf$y / sqrt(c)
    # return a list with the data
    list(js = js, dJdt = dJdt, wf = s$wf, u2j = s$u2j, name = paste0(kernelName, '.', n))
  }

  findReggeonData <- cache(findReggeonDataFun, kernelName)

  # find the complete kernel
  findKernel <- function(.t = 0, ...) {
    # use parallelization if possible
    lapplyType <- if(Sys.info()['sysname'] == 'Linux') mclapply else lapply
    lapplyType(1:numReg, findReggeonDataFun, .t, ...)
  }

  k <- list(findKernel = findKernel,
            findReggeonData = findReggeonData,
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
