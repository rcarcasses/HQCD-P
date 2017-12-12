# data is a list of the form (F2, Q2, x, err)
#' @export
F2model <- function(potPars = 'bcde',
                    numReg = 3,
                    data = NULL,
                    addAlfonsoTerm = TRUE,
                    considerNonMinimalCoupling = TRUE,
                    forceSPjsWeight = 10^5,
                    ihqcdPars = ''
                    ) {
  if(is.null(data))
    data <- getHERAF2()

  # this variable contains all the kernels used in the model
  modelKernels  <- c()
  colNames      <- c()
  # the parameters to be optimized
  # we should keep in mind that there may be shared
  # parameters between different kernels, like ls
  # so by using its name as identifier we effectively
  # acomplish that
  parameters    <- c()

  # add a new kernel to the fitting F2  model
  addKernel <- function(...) {
    kernelData <- list(...)
    # create a f2 object with the configuration desired
    f2 <- F2(numReg = kernelData$numReg,
             data = data,
             kernelName = kernelData$prefix,
             considerNonMinimalCoupling = considerNonMinimalCoupling)
    # configure the f2 object with the custom potential passed
    f2$setNewPotential(kernelData$potential, kernelData$optimPars)
    # add the created object to the kernelData one
    kernelData <- c(kernelData, f2 = f2, list())
    # add this kernel to the list of kernels
    modelKernels[[kernelData$prefix]] <<- kernelData
    # update the parameters vector
    mapply(function(n, v) parameters[n] <<- v, names(kernelData$optimPars), kernelData$optimPars)

    # and update the list of column names
    colNames <<- lapply(modelKernels, function(k) {
                   # create a vector with the names = prefix + index
                   ns <- lapply(0:(k$numReg - 1), function(i) paste0(k$prefix, i))
                   if(considerNonMinimalCoupling) {
                     # add the NMC names
                     ns <- c(ns, lapply(0:(k$numReg - 1), function(i) paste0(k$prefix, 'NMC', i)))
                   }
                   # return the names vector created
                   ns
                 })
  }

  # prints the information of the current state of the model
  info <- function() {
    cat('\n\t##############################')
    cat('\n\t# Overview of the model')
    cat('\t- Number of kernels', length(modelKernels), '\n')
    cat('\t- Reggeons by kernel', paste0(unlist(lapply(modelKernels, function(k) k$numReg)), collapse = ', '), '\n')
    cat('\t- Column names', getColNames(), '\n')
    cat('\t- Parameters to optimize are \n')
    print(parameters)
    cat('\n\t##############################')
    cat('\n')
  }

  # returns a vector of strings with the current column names
  getColNames <- function() unlist(lapply(colNames, function(col) unlist(col)))

  ihqcdPars <- strsplit(ihqcdPars, split = ',')[[1]] # convert the string to a vector
  lastEval  <- NULL
  bestEval  <- NULL


  rssFunBase <- function(pars) {
    dfCmdString <- paste0('data.frame(',
                            paste0(colNames, collapse = ', '),
                          ')')
    df <- eval(parse(text = dfCmdString))

    # now get the best linear fit for the coefficients
    fsString <- paste0(colNames, collapse = ' + ')
    # this returns a list with the values of the best coefficients found
    lmString <- paste0('lm(data$F2~', fsString,' - 1, weights = 1/data$err^2)$coefficients')
    lmString
    # coeffs <- eval(parse(text = lmString))
  }

  # compute the Residual Sum of Squares for the
  # given configuration of kernels and the parameters
  # passed. The parameters here are the one at the
  # potentials, couplings are deduced from a linear fit.
  rssFunBase <- function(pars) {
    # dynamic parameters will be sent
    # the name of the parameters are lost, here we fix that
    names(pars) <- names(parameters)

    # now we compute each one of the fns for each kernel
    df <- NULL
    lapply(modelKernels, function(k) {
      # get the parameters associated with this kernel
      indices <- names(k$optimPars)
      kPars <- pars[indices]
      # get the right function to call according whether we are considering
      # the non minimal coupling case or not
      f <- ifelse(considerNonMinimalCoupling, k$f2.getAllFns, k$f2.getFns)
      # get the dataframe with the fns
      kdf <- do.call(f, as.list(kPars))
      # add the computed dataframe columns to the full one
      if(is.null(df))
        df <<- kdf
      else
        df <<- cbind(df, kdf)
    })
    # now with the full dataframe of fns we proceed to find the
    # best linear fit for the coefficients multiplying this
    colnames(df) <- getColNames()

    # now get the best linear fit for the coefficients
    fsString <- paste0(unlist(lapply(colnames(df), function(n) paste0('df$', n ))), collapse = ' + ')
    # this returns a list with the values of the best coefficients found
    lmString <- paste0('lm(data$F2 ~ ', fsString,' - 1, weights = 1/data$err^2)')
    lm <- eval(parse(text = lmString))
    coeff <- lm$coefficients
    print(coeff)
    s    <- sum((resid(lm)/data$err)^2)
    chi2 <- s / (length(data$F2) - length(parameters))
    # cat('\r                                                                => chi2 =', chi2)

    # always remember the last evaluation
    lastEval <<- list(pars = pars, coeff = coeff, val = s, chi2 = chi2)
    # and also keep track of the best result found
    if(!is.na(s)) {
      if(is.null(bestEval)) {
        bestEval <<- list(pars = pars, coeff = coeff, val = s, chi2 = chi2)
        cat('\nstarting chi2 =', chi2)
        rredis::redisSAdd(potPars, bestEval)
      } else if(bestEval$val > s ) {
        bestEval <<- list(pars = pars, coeff = coeff, val = s, chi2 = chi2)
        flog.debug(paste('\t *chi2* =', chi2))
        flog.debug(paste('\t  coeff =', coeff, ''))
        rredis::redisSAdd(potPars, bestEval)
      }
    }

    # we want to force the soft pomeron intercept to 1.09, here is where we do it
    # we will take the second intercept from the first kernel added and
    # will add a penalization term proportional to the difference to the desired value
    spJs <- modelKernels[[1]]$f2.getJs()[2]
    flog.debug(paste('intercepts ', modelKernels[[1]]$f2.getJs(), ' spJs =', spJs))
    s <- s + forceSPjsWeight * (spJs - 1.09)^2
    s
  }

  rssFun <- function(pars) {
    rssFunBase(pars)
  }

  chain <- function(method) {
    flog.debug(paste('[F2model] Using last evaluation parameters\n', bestEval$pars))
    startOpt <<- bestEval$pars
    bestEval <<- NULL
    optimFun(method)
  }

  # call the optimization function as long as it gives a better result
  # if the difference between two successive attempts is less than tol
  # then it stops
  fit <- function(method = NULL, tol = 1e-3) {
    op <- NULL
    i <- 1
    callFun <- function(fn) {
      tic()
      op <<- fn(method)
      exectime <- toc()
      exectime <- exectime$toc - exectime$tic
      flog.debug(paste('[F2model]  -::- Iteration', i, ' completed in', round(exectime / 60), ' min, new chi2', bestEval$chi2, ' ---'))
      i <<- i + 1
    }

    if(is.null(bestEval))
      callFun(optimFun)

    if(is.null(op))
      return(NULL)

    while(TRUE) {
      lastBestChi2 <- bestEval$chi2
      callFun(chain)
      newBestChi2 <- bestEval$chi2
      if(abs(newBestChi2 - lastBestChi2) < tol) # the new iteration wasn't better than the one before
        break

      lastBestChi2 <- newBestChi2
    }
    op
  }

  optimFunRaw <- function(method = NULL) {
    if(is.null(method))
      method <- "Nelder-Mead"

    flog.debug(paste('[F2model] Using method', method))

    if(method == 'hjk')
      op <- hjk(startOpt, fn = rssFun)
    else if(method == 'sbplx')
      op <- sbplx(startOpt, fn = rssFun)
    else if(method == 'newuoa')
      op <- newuoa(startOpt, fn = rssFun)
    else if(method == 'nmk')
      op <- nmk(startOpt, fn = rssFun)
    else if(method == 'cmaes')
      op <- cma_es(startOpt, fn = rssFun)
    else
      # go to default optim implementation
      op <- optim(parameters,
          fn = rssFun,
          method = method,
          hessian = FALSE#,
          #upper = upperOpt,
          #lower = lowerOpt,
          #control = list(trace = 10, REPORT = 500, maxit = 1e6)
          )
  }

  optimFun <- optimFunRaw

  alreadyComputed <- function() {
    cat('Coming soon!\n')
  }

  plotF2 <- function() {
    flog.error(paste('[F2model] Update to use the new kernels configuration'))
    cat('Using', bestEval$pars, bestEval$coeff,'\n')
    do.call(f2$plot, as.list(c(parsSet, bestEval$pars, bestEval$coeff)))
  }

  plotSpectrum <- function() {
    flog.error('[F2model] Update to use the new kernels configuration')
    flog.debug(paste('Using', bestEval$pars, '\n'))
    do.call(f2$plotSpectrum, as.list(c(parsSet, bestEval$pars)))
  }

  plotRegge <- function(glueballs = TRUE, mesons = TRUE) {
    flog.error('[F2model] Update to use the new kernels configuration')
    flog.info(paste('Using', bestEval$pars))
    do.call(f2$plotRegge, as.list(c(parsSet, bestEval$pars)))

    if(glueballs)
      f2$plotGlueballMasses()

    if(mesons)
      f2$plotMesonsMasses()
  }

  plotEffectiveExponent <- function() {
    flog.error('[F2model] Update to use the new kernels configuration')
    do.call(f2$plotEffectiveExponent, as.list(c(bestEval$pars, bestEval$coeff)))
  }

  a <- list(optimFun = optimFun,
            chain = chain,
            fit = fit,
            info = info,
            addKernel = addKernel,
            plotF2 = plotF2,
            plotRegge = plotRegge,
            plotSpectrum = plotSpectrum,
            plotEffectiveExponent = plotEffectiveExponent,
            getColNames = function() colNames,
            getData = function() data,
            getLastEval = function() lastEval,
            getBestEval = function() bestEval,
            getParsEval = function() parsEval,
            getParsSet = function() parsSet,
            getPotPars = function() potPars)

  class(a) <- append(class(a), 'attempt')
  a
}
