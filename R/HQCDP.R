#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function(alpha = 0) {
  h <- list(processes = list(), kernels = list())
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
  # add the constraint for the intercept of the soft pomeron
  # the value of this attribute will be used as weight while fitting
  attr(h, 'addSPconstraint') <- 1e6
  attr(h, 'gtOrder')         <- 1
  attr(h, 'alpha')           <- alpha
  h
}

#' @export
getNeededTVals.HQCDP <- function(p) sort(unique(c(0, unlist(lapply(p$processes, getNeededTVals)))))
#' @export
enlargeData.HQCDP <- function(p) lapply(p$processes, enlargeData)

#' @export
addKernel <- function(x, ...) UseMethod('addKernel')
#' @export
addKernel.default <- function(x, ...) 'calling addKernel in the wrong object'
#' @export
addKernel.HQCDP <- function(h, ...) {
  k <- kernelUnit(...)
  h$kernels <- append(h$kernels, list(k))
  # return the modified HQCD object
  h
}
#' @export
addProcessObservable <- function(x, ...) UseMethod('addProcessObservable')
#' @export
addProcessObservable.default <- function(x, ...) 'calling addProcessObservable in the wrong object'
#' @export
addProcessObservable.HQCDP <- function(h, x) {
  # insert the non minimal coupling attribute in the child processes
  attr(x, 'alpha') <- attr(h, 'alpha')
  if(is.element('ProcessObservable', class(x)))
    h$processes <- append(h$processes, list(x))
  else
    flog.error(paste('Adding an object of unknown classes:', class(x)))
  # return the modified HQCD object
  h
}

#' @export
print.HQCDP <- function(x) {
  cat('HQCD processes:', length(x$processes), ' kernels:', length(x$kernels), '\n')
}

#' This computes the sum of the rss for all the processes configured.
#' This is the function to minimize in the optimization procedure.
#' @export
rss.HQCDP <- function(x, pars = NULL, allGs = NULL, startGs = NULL) {
	#cat('\n allGs', allGs, ' startGs', startGs, '\n')
  # warn the user that there are configurations
  # that are still required
  if(length(x$processes) == 0 || length(x$kernels) == 0) {
    cat('Please complete your HQCDP configuration first\n')
    return()
  }
  # first we need to compute the spectra of the kernels
	# this is a parallelized call
  spectra <- getSpectra(x, pars)
  # now we find all the fns for each one of the processes
  allProcFns <- lapply(x$processes, getFns, spectra = spectra)
	gs <- NULL
	# if allGs is set then just evaluate the function
	val <- if(!is.null(allGs)) {
		gs <- allGs
		evalRSSInGs(x, allProcFns, allGs)
	}
	else {
		# else find the best gs by going through the gs submanifold optimization
		bestGs <- getBestGs(x, allProcFns, startGs)
		gs     <- unlist(bestGs$gs)
		bestGs$value
	}
  # add SP constraint here
  valWeighted <- 0
  if(!is.null(attr(x, 'addSPconstraint'))) {
    spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
    jsSP <- spectraForTZero[[1]][[2]]$js
    # cat('SP intercept', jsSP, '\n')
    valWeighted <- val + attr(x, 'addSPconstraint') * (jsSP - 1.09)^2
  }
  # return the result in the appropiated format
  if(is.null(attr(x, 'complete')))
    val
  else
    list(val = val, valWeighted = valWeighted, gs = gs)
}

#' Given an ProcessObservable or HQCDP object and its associated fns already computed
#' find the best coefficients, gs, such that the sum of difference squared
#' weigthed by the inverse of the error squared (rss) between the prediction
#' and the experimental values of the process is minimized.
#' @export
getBestGs <- function(x, ...) UseMethod('getBestGs')
#' @export
getBestGs.default <- function(x) paste('getBestGs has to be implemented for this object with classes', class(x))

#' Given a pre computed processes fns find the best values for the gs
#' @return an optim object
#' @export
getBestGs.HQCDP <- function(x, allProcFns, startGs) {
  # first we need to define an function depending only of the gs
  # to be optimized
	fn <- function(parAllGs) evalRSSInGs(x, allProcFns, parAllGs)
	numGs <- (attr(x, 'gtOrder') + 1) * sum(unlist(lapply(x$kernels, '[[', 'numReg')))
	if(attr(x, 'alpha') != 0)
	  numGs <- 2 * numGs
  # if is the fist time just put something there
  if(is.null(startGs))
    startGs <- rep(1, len = numGs)
	# define the global gradient function
  grad <- function(allGs) {
    gs <- gs.as.data.frame(x, allGs)
    rowSums(as.data.frame(mapply(function(proc, procFns) {
      gradRSSGs(proc, fns = procFns, gs = gs)
    }, x$processes, allProcFns)))
  }
  op <- optim(startGs,
						 	#gr = grad,
						 	fn = fn, method = 'BFGS', hessian = FALSE, control = list(maxit = 1000))
  # store the best gs found so they can be used as a starting point of the next call
	flog.debug(paste('bestGs  =', do.call(paste, as.list(format(op$par, digits = 4))), ' in', op$counts[1], ' steps'))
  # add the property gs with the dataframe format
  op$gs <- gs.as.data.frame(x, op$par)
  op
}

gs.as.data.frame <- function(x, allGs) {
  as.data.frame(matrix(allGs, ncol = attr(x, 'gtOrder') + 1))
}

# given all processes fns and a vector of gs compute the rss sum
evalRSSInGs <- function(x, allProcFns,  allGs) {
	# build a dataframe of gs, which is what the predict function of the ProcessObservable
	# is expecting, from the allGs passed
	# pay attention to the number of columns: is related to the g(t) order expansion
	gs <- gs.as.data.frame(x, allGs)
	# find the rss for each process for the given values of gs and fns
	sum(unlist(mapply(function(proc, procFns) {
		# remove some possible NA from the fns and put some large number
		if(any(is.na(procFns))) {
			flog.warn('There were NA values for this evaluation')
			1e30
		}
		else
			rss(proc, fns = procFns, gs = gs)
	}, x$processes, allProcFns)))
}

#' @export
fit <- function(x, ...) UseMethod('fit')
#' @export
fit.HQCDP <- function(x, allPars = NULL, initGs = NULL, method = 'Nelder-Mead') {
	flog.debug('Using method %s', method)
  DoF <- getDoF(x)
	flog.debug(' --- DoF %s', DoF)
  # here we declare that we want the full output of the rss function for instance
  attr(x, 'complete') <- TRUE
  # get all the parameters to fit if required
  if(is.null(allPars))
    allPars <- getKernelPars(x)
  # reset the bestEvalEnv
  resetEvalTracker(allPars)
  # the function to internally called by optim
  fn <- function(pars) {
    names(pars) <- names(allPars)
		# fn arguments, allGs set only the first time, then use last gs as starting point
		# this is just to allow reproducibility of known results at the first ru
    lastEval    <- get('lastEval', envir = bestEvalEnv)
    completeVal <- rss(x, pars, allGs = if(is.null(lastEval$gs)) initGs else NULL, startGs = lastEval$gs)
		val         <- completeVal$val
    valWeighted <- completeVal$valWeighted
    gs          <- completeVal$gs
    chi2        <- val / DoF
    flog.debug('          pars = %s', do.call(paste, as.list(round(pars, 6))))
		flog.debug('          chi2 = %s', round(chi2, 3))
    # store the partial results in the best eval tracker
    saveStep(chi2, val, pars, gs)
		# optimize the log better
    log(valWeighted)
  }

  i <- 1
  while(TRUE) {
    tic()
    bestEval <- get('bestEval', envir = bestEvalEnv)
    lastBestChi2 <- bestEval$chi2
    startPars <- bestEval$pars
    op <- optim(startPars,
								fn = fn,
								hessian = FALSE, method = method)
    bestEval <- get('bestEval', envir = bestEvalEnv)
    newBestChi2 <- bestEval$chi2
    if(abs(newBestChi2 - lastBestChi2) < 1e-3) # the new iteration wasn't better than the one before
      break

    lastBestChi2 <- newBestChi2
    exectime <- toc()
    flog.debug(paste('[HQCDP]  -::- Iteration', i, ' completed in', seconds_to_period(round(exectime$toc - exectime$tic)), ' min, new chi2', bestEval$chi2))
    i <<- i + 1
  }
  op
}

getKernelPars <- function(x) {
  allPars <- unlist(lapply(x$kernels, function(k) k$optimPars))
  keep    <- !duplicated(names(allPars))
  allPars[keep]
}

getDoF <- function(x) {
  # get the experimental points per process
  expPoints <- sum(unlist(lapply(x$processes, function(p) length(expVal(p)))))
  # get the amount of fitting parameters
  # gDeg * number of Reggeons since g(t) = g0 + g1 * t +...
	numGs <- (attr(x, 'gtOrder') + 1) * sum(unlist(lapply(x$kernels, '[[', 'numReg')))
  fitParams <- numGs + length(getKernelPars(x))
  expPoints - fitParams
}

# store the number of cores for this computation
cores <- if(Sys.getenv('USE_CORES') == '') {
	# Calculate the number of cores, left one for the system
	detectCores() - 1
} else {
	# use the value in the USE_CORES system environment variable
	# as the amount of desired cores
	as.integer(Sys.getenv('USE_CORES'))
}

flog.debug('Using %s cores', cores)
#' @export
getSpectra <- function(x, ...) UseMethod('getSpectra')
#' @export
getSpectra.default <- function(x) 'getSpectra not defined in current object'
#' Get the spectra of a configured HQCDP object
#' @param x a HQCDP object, already configured
#' @param ... additional parameters which will be passed
#' down while computing the spectrum
#' @export
getSpectra.HQCDP <- function(x, pars = NULL, ts = NULL) {
  if(is.null(ts))
	  ts <- getNeededTVals(x)

	numRegs <- unlist(lapply(x$kernels, function(k) k$numReg))
	# unwrap the computation of each one of the Reggeon data for each one of the kernels
	unwrappedFunCalls <- Reduce(function(acct, t) {
  	as.list(c(acct,
  	  Reduce(function(acc, k) {
    	  as.list(c(acc,
    	    lapply(1:k$numReg, function(n) {
    	      kArgs <- c(list(n = n, .t = t), pars)
    	      # return the function to be called
    	      function(i) {
              # we need to initialize the computation on each node
              init()
              do.call(k$findReggeonData, kArgs)
            }
    	    }))
    	  )
  	  }, x$kernels, init = list()))
  	)
	}, ts, init = list())

  # now compute every kernel for each value of t and the parameters passed
  # this is the most expensive part of the computation and its parallelizable
  # so we use parLapply. Notice that we are already using mclapply in each of the
  # kernel computations so this is the most we can go while nesting parallel computations
  convertRawSpectra(
    if(Sys.info()['sysname'] == 'Linux')
      mclapply(unwrappedFunCalls, function(f) f(), mc.cores = cores)
    else
      lapply(unwrappedFunCalls, function(f) f())
  , numRegs, ts)
}

#' @export
convertRawSpectra <- function(rawSpectra, numRegs, ts) {
  tBlockSize <- length(rawSpectra) / length(ts)
  lapply(1:length(ts), function(i) {
    t     <- ts[i]
    stRaw <- rawSpectra[(1:tBlockSize) + (i - 1) * tBlockSize]
    # mark the slots, output is like
    # [1] 1 1 2 2 2 3
    slots <- unlist(lapply(1:length(numRegs), function(n) rep_len(n, numRegs[n])))
    st <- lapply(1:length(numRegs), function(n) {
      positions <- n == slots
      stRaw[positions]
    })
    list(t = t, spectra = st)
  })
}

#' @export
plot.HQCDP <- function(x, pars = NULL, gs = NULL) {
  # we need to compute the spectra for an enlarged set of t values
  # convert gs to a data frame, if required
  if(!is.data.frame(gs))
    gs <- gs.as.data.frame(x, gs)
  # get the plot points
  plotPoints <- enlargeData(x)
  # compute the spectrum, now with the particular t values needed
  ts <- sort(unique(unlist(lapply(plotPoints, function(df) if(!is.null(df$t)) df$t))))
  pb <- txtProgressBar(min = 0, max = 100, initial = 1, style = 3)
  spectra <- getSpectra(x, pars, ts)
  setTxtProgressBar(pb, 30)
  # get the fns for the plot points
  i <- 0
  allProcFns <- mapply(function(proc, points) {
                        setTxtProgressBar(pb, 30 + 30 * i / length(x$processes))
                        i <<- i + 1
                        list(getFns(proc, spectra = spectra, points = points))
                      }, x$processes, plotPoints)
  # find the predictions for the plot points for each one of the processes
  i <- 0
  predPlotPoints <- mapply(function(proc, procFns, procPlotPoints) {
    pred <- predict(proc, points = procPlotPoints, fns = procFns, gs = gs)
    setTxtProgressBar(pb, 60 + 40 * i / length(x$processes))
    i <<- i + 1
    #cat('predictions found', unlist(pred), '\n')
		list(cbind(procPlotPoints, data.frame(predicted = pred)))
	}, x$processes, allProcFns, plotPoints)
  # call the plot function on each on
  mapply(plot, x$processes, predPlotPoints)
  # return invisibly the computed points, useful to speed up future plots
  invisible(predPlotPoints)
}
