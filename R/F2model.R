# data is a list of the form (F2, Q2, x, err)
#' @export
F2model <- function(potPars = 'bcde',
                    numReg = 3,
                    data = NULL,
                    addAlfonsoTerm = TRUE,
                    considerNonMinimalCoupling = TRUE,
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
    colNames <<- lapply(modelKernels, function(k)
                      lapply(0:(k$numReg - 1), function(i)
                                                paste0(k$prefix, i)))
  }

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
      # get the dataframe with the fns
      kdf <- do.call(k$f2.getFns, as.list(kPars))
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
    # multiply each colum of the data frame by the correspondent coefficient
    # we need to transpose the data frame twice to achieve this
    #df <- t(t(df) * coeff)
    # to get the estimated F2 we need to sum all the rows
    #F2df <- rowSums(df)
    #chi2vec <- unlist(mapply(function(F2, F2pred, err)
    #                          ((F2pred - F2) / err)^2, data$F2, F2df, data$err))
    #s    <- sum(chi2vec)
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
        cat('\n\t *chi2* =', chi2)
        cat('\n\t  coeff =', coeff, '')
        rredis::redisSAdd(potPars, bestEval)
      }
    }

    # we want to force the soft pomeron intercept to 1.09, here is the right place to do it
    #s <- s + 1000000 * rss$size * (rss$Js[[2]] - 1.09)^2
    cat('\nTODO: force soft pomeron intercept\n')
    s
  }

  rssFun <- function(pars) {
    rssFunBase(pars)
  }

  chain <- function(method) {
    cat('\nUsing last evaluation parameters\n', bestEval$pars, '\n')
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
      cat('\n --- Iteration', i, ' completed in', round(exectime / 60), ' min, new chi2', bestEval$chi2, ' ---\n')
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

    cat('\n Using method', method, '\n')

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
    cat('[WARN] Update to use the new kernels configuration\n')
    cat('Using', bestEval$pars, bestEval$coeff,'\n')
    do.call(f2$plot, as.list(c(parsSet, bestEval$pars, bestEval$coeff)))
  }

  plotSpectrum <- function() {
    cat('[WARN] Update to use the new kernels configuration\n')
    cat('Using', bestEval$pars, '\n')
    do.call(f2$plotSpectrum, as.list(c(parsSet, bestEval$pars)))
  }

  plotRegge <- function(glueballs = TRUE, mesons = TRUE) {
    cat('[WARN] Update to use the new kernels configuration\n')
    cat('Using', bestEval$pars, '\n')
    do.call(f2$plotRegge, as.list(c(parsSet, bestEval$pars)))

    if(glueballs)
      f2$plotGlueballMasses()

    if(mesons)
      f2$plotMesonsMasses()
  }

  plotEffectiveExponent <- function() {
    cat('[WARN] Update to use the new kernels configuration\n')
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
