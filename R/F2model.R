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
  modelKernels <- c()
  colNames     <- c()
  # the parameters to be optimized
  # we should keep in mind that there may be shared
  # parameters between different kernels, like ls
  # so by using its name as identifier we effectively
  # acomplish that
  parameters   <- c()
  bestEval     <- NULL
  startOpt     <- NULL

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
    kernelData <- c(kernelData, list(f2 = f2))
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

  # get the Fns, with x factor included, for each one of the
  # kernels configured in this run. The order is the same as
  # the one used for adding the kernels and can be tracked
  # through the column names of the returned dataframe
  getFns <- function(pars) {
    df <- NULL
    lapply(modelKernels, function(mk) {
      # get the parameters associated with the potential of this kernel
      indices <- names(mk$optimPars)
      mkPars <- as.list(pars[indices])
      # if there is a data object, then clone it as well
      # this is used to compute Fns for custom values of Q2 and x
      if(is.list(pars) && !is.null(pars$data))
        mkPars$data <- pars$data

      # get the right function to call according whether we are considering
      # the non minimal coupling case or not
      f <- if(considerNonMinimalCoupling) mk$f2$getAllFns else mk$f2$getFns
      # get the dataframe with the fns
      mkdf <- do.call(f, mkPars)
      # add the computed dataframe columns to the full one
      if(is.null(df))
        df <<-mkdf
      else
        df <<- cbind(df, mkdf)
    })
    # put the right column names
    colnames(df) <- getColNames()
    df
  }

  # given a dataframe with the Fns at its column, find the
  # best linear combination of Fns to explain experimental data
  getLinearFnsFit <- function(df) {
    fsStr <- paste0(unlist(lapply(colnames(df), function(n) paste0('df$', n ))), collapse = ' + ')
    # this returns a list with the values of the best coefficients found
    lmStr <- paste0('lm(data$F2 ~ ', fsStr,' - 1, weights = 1/data$err^2)')
    eval(parse(text = lmStr))
  }

  runOnce <- function(pars) {
    # the name of the parameters are lost, here we fix that
    names(pars) <- names(parameters)
    # now we compute each one of the fns for each kernel
    df <- getFns(pars)
    # now with the full dataframe of fns we proceed to find the
    # best linear fit for the coefficients multiplying this
    lm    <- getLinearFnsFit(df)
    coeff <- lm$coefficients
    rss   <- sum((resid(lm)/data$err)^2)
    chi2  <- rss / (length(data$F2) - length(parameters))
    # evaluation result
    value <- list(pars = pars,
                  coeff = coeff,
                  rss = rss,
                  chi2 = round(chi2, digits = 4),
                  df = df,
                  jsk1 = round(modelKernels[[1]]$f2$getJs(), digits = 3))
    # and keep track of the best result found
    if(!is.na(rss)) {
      checkBestEval(value)
    }
    value
  }

  # this function's goal is to allow the user to make a single run
  # with the parameters specified and therefore to setup the best
  # evaluation without making the full optimization.
  singleRun <- function(pars = NULL) {
    if(is.null(pars))
      pars <- parameters

    # do not return anything
    v <- runOnce(pars)
  }

  # compute the Residual Sum of Squares for the
  # given configuration of kernels and the parameters
  # passed. The parameters here are the one at the
  # potentials, couplings are deduced from a linear fit.
  rssFun <- function(pars) {
    runStep     <- runOnce(pars)
    rssWeighted <- runStep$rss
    # we want to force the soft pomeron intercept to 1.09, here is where we do it
    # we will take the second intercept from the first kernel added and
    # will add a penalization term proportional to the difference to the desired value
    spJs        <- runStep$jsk1[2]
    rssWeighted <- rssWeighted + forceSPjsWeight * (spJs - 1.09)^2
    rssWeighted
  }

  # record the evaluation with the best chi2
  checkBestEval <- function(value) {
    if(is.null(bestEval)) {
      # first time
      bestEval <<- value
      flog.info('starting chi2 = %s', bestEval$chi2)
      flog.debug(paste('\t\t jsk1 =', do.call(paste, as.list(bestEval$jsk1))))
    } else if(bestEval$rss > value$rss ) {
      bestEval <<- value
      flog.info(paste('\t *chi2* =', bestEval$chi2, 'pars', do.call(paste, as.list(round(bestEval$pars, digits = 3)))))
      flog.debug(paste('\t\t jsk1 =', do.call(paste, as.list(bestEval$jsk1))))
    }
  }

  # predicts the values of the model for a given set of Q2 and x
  # the predData should be a dataframe with the same structure as the
  # one returned by getHERAF2() function, currently
  # the Q2 are restricted to the experimental values only
  predictF2 <- function(predData = NULL) {
    # if no data is passed, then use the experimental one
    if(is.null(predData))
      predData <- data

    # force the user to fit the model first
    if(is.null(bestEval)) {
      cat('Train the model first, try fit() or singleRun() functions \n')
      return()
    }
    flog.debug('[F2model] Predicting values of F2 for best parameters found %s', dumpList(bestEval$pars))
    # get the Fns for the best parameters found
    # specify the data to produce the Fns
    df <- getFns(c(bestEval$pars, list(data = predData)))
    # predict the values from the Fns and the coefficients
    apply(df, 1, function(fn) sum(fn * bestEval$coeff))
  }

  chain <- function(method) {
    flog.debug(paste('[F2model] Using last evaluation parameters', do.call(paste, as.list(bestEval$pars))))
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
      flog.debug(paste('[F2model]  -::- Iteration', i, ' completed in', seconds_to_period(round(exectime / 60)), ' min, new chi2', bestEval$chi2))
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
    if(is.null(startOpt))
      startOpt <<- parameters

    if(is.null(method))
      method <- "Nelder-Mead"

    flog.debug(paste('[F2model] Using method', method, ' starting at ', do.call(paste, as.list(startOpt))))

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
      op <- optim(startOpt,
                  fn = rssFun,
                  method = method,
                  hessian = FALSE)
    op
  }

  optimFun <- optimFunRaw

  alreadyComputed <- function() {
    cat('Coming soon!\n')
  }

  # plot the F2 experimental data together
  # with the prediction of this model.
  # extra parameters will be passed to the plot function
  plotF2 <- function(...) {
    flog.debug(paste('[F2model] Ploting F2 for the last model parameters'))
    plotHERAF2(...)
  }

  plotHERAF2 <- function(...) {
    if(is.null(bestEval)) {
      cat('Train the model first, try fit() or singleRun() functions \n')
      return()
    }
    # plot the experimental points
    xLimMax <- 2 * max(data$x)
    plotPars <- c(list(x = data$x, y = data$F2, type = "n", log = 'x',
                       xlim = c(1e-6, xLimMax), ylim = c(0.6 * min(data$F2), 1.2 * max(data$F2)),
                       xlab = expression("x"), ylab = expression("F"[2])), xaxt = "n", list(...))
    do.call(plot, plotPars)
    hTicks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
    vTicks <- seq(0.5, 2, len = 4)
    abline(v = hTicks, h = vTicks, col = 'gray90', lty = 3)
    axis(1, at = hTicks)
    axis(2, at = vTicks)
    # create a dataframe from hera data
    df <- data.frame(F2 = data$F2, Q2 = data$Q2, x = data$x, err = data$err)
    # add some alternating colors
    df <- addAlternatingColToDataByQ2(df, colName = 'color', col = c('blue', 'red', 'green'))
    # add some alternating shapes for the scattered plot
    df <- addAlternatingColToDataByQ2(df, colName = 'pch', col = 0:6)
    # do the plot
    Q2s <- unique(df$Q2)
    # progress bar
    pb <- txtProgressBar(min = 0, max = length(Q2s), style = 3)
    lapply(Q2s, function(Q2) {
      Q2data <- df[df$Q2 == Q2,]
      # plot the experimental points
  	  lines(Q2data$x, Q2data$F2, type = "p", col = Q2data$color, cex = 0.5, pch = Q2data$pch)
  	  # create a vector of x that covers the experimental points
  	  x <- c(0.5 * min(Q2data$x))
  	  # put 10 theoretical points between two experimental
  	  lapply(Q2data$x, function(xVal) x <<- c(x, seq(x[length(x)], xVal, len = 10)))
  	  # and add the end point
  	  x <- c(x, 2 * max(Q2data$x))
  	  # now compute the predicted values of F2 for the above values of x and the current Q2
  	  predF2 <- predictF2(predData = list(Q2 = rep_len(Q2, length(x)), x = x))
  		# plot the dots
      lines(x, predF2, type = 'l', col = Q2data$color)
  		# plot the error bars
  		arrows(Q2data$x, Q2data$F2 - Q2data$err, Q2data$x, Q2data$F2 + Q2data$err,
  		       length = 0.02, angle = 90, code = 3, col = Q2data$color)

  		# put the Q2 value at the end of the lined
  		textPosX <- 0.35 * min(Q2data$x)
  		#textPosX <- 10^(0.65 * (log(x, base = 10) + 5))
  		#cat('textX', textPosX,'\n')
  		textPosY <- splinefun(x, predF2)(textPosX)
      boxed.labels(textPosX, textPosY, labels = paste(Q2), col = Q2data$color, cex = 0.7, xpad = 1, ypad = 1, border = FALSE, bg = 'white')
      # update the progress bar
      setTxtProgressBar(pb, match(Q2, Q2s))
    })

    # empty return
    return()
  }

  plotSpectrum <- function() {
    lapply(modelKernels, function(mk) {
      # call plotSpectrum on each kernel object with the best
      # parameters found
      do.call(mk$f2$plotSpectrum, as.list(bestEval$pars))
    })
    flog.info('[F2model] Spectrum for each kernel plotted')
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
            singleRun = singleRun,
            predict = predictF2,
            info = info,
            addKernel = addKernel,
            plot = plotF2,
            plotRegge = plotRegge,
            plotSpectrum = plotSpectrum,
            plotEffectiveExponent = plotEffectiveExponent,
            getColNames = function() colNames,
            getData = function() data,
            getBestEval = function() bestEval,
            getPotPars = function() potPars)

  class(a) <- append(class(a), 'attempt')
  a
}
