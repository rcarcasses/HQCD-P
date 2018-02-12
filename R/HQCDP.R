#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function() {
  h <- list(processes = list(), kernels = list(), cluster = NULL)
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
  # add the constraint for the intercept of the soft pomeron
  # the value of this attribute will be used as weight while fitting
  attr(h, 'addSPconstraint') <- 1e6
  h
}

#' @export
setUpCluster <- function(x, ...) UseMethod('setUpCluster')
setUpCluster.default <- function(x, ...) 'setUpCluster not implemented for this object'
#' @export
setUpCluster <- function(x, ncores) {
  x$cluster <- makeCluster(ncores, outfile = 'cluster.log')
  # we need to put the ihqcd information available in the cluster
  clusterExport(x$cluster, names(solve(iHQCD())))
  x
}

#' @export
destroyCluster <- function(x, ...) UseMethod('destroyCluster')
destroyCluster.default <- function(x, ...) 'destroyCluster not implemented for this object'
#' @export
destroyCluster <- function(x, ncores) {
  if(!is.null(x$cluster)) {
    flog.debug('Stopping cluster\n')
    parallel::stopCluster(x$cluster)
  }
  x$cluster <- NULL
  x
}

#' @export
getNeededTVals.HQCDP <- function(p) unique(unlist(lapply(p$processes, getNeededTVals)))

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
rss.HQCDP <- function(x, ...) {
  # warn the user that there are configurations
  # that are still required
  if(length(x$processes) == 0 || length(x$kernels) == 0) {
    cat('Please complete your HQCDP configuration first\n')
    return()
  }
  # first we need to compute the spectra of the kernels
	# this is a parallelized call
  spectra <- getSpectra(x, ...)
  # now we find all the fns for each one of the processes
  allProcFns <- lapply(x$processes, getFns, spectra = spectra)
  # now we proceed to the gs submanifold optimization
  bestGs <- getBestGs(x, allProcFns)
  val <- bestGs$value
  #val <- sum(unlist(lapply(x$processes, rss, spectra = spectra)))
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
    list(val = val, valWeighted = valWeighted, gs = bestGs$gs)
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
getBestGs.HQCDP <- function(x, allProcFns) {
  # first we need to define an function depending only of the gs
  # to be optimized
  fn <- function(allGs) {
    # build a dataframe of gs, which is what the predict function of the ProcessObservable
    # is expecting, from the allGs passed
    # pay attention to the number of columns
    gs <- gs.as.data.frame(allGs)
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
  # now use some good starting point for this optimization
  startGs <- get('startGs', envir = bestEvalEnv)
	flog.debug(paste('getBestGs: startGs =', do.call(paste, as.list(format(startGs, digits = 4)))))
  # if is the fist time just put something there
  if(is.null(startGs))
    startGs <- rep(1, len = 2 * sum(unlist(lapply(x$kernels, '[[', 'numReg'))))

  grad <- function(allGs) {
    gs <- gs.as.data.frame(allGs)
    rowSums(as.data.frame(mapply(function(proc, procFns) {
      gradRSSGs(proc, fns = procFns, gs = gs)
    }, x$processes, allProcFns)))
  }
  # TODO: chain optimizations (?)
  op <- optim(startGs, gr = grad, fn = fn, method = 'BFGS', hessian = FALSE)
  # store the best gs found so they can be used as a starting point of the next call
  assign('startGs', op$par, envir = bestEvalEnv)
  # add the property gs with the dataframe format
  op$gs <- gs.as.data.frame(op$par)
  op
}

gs.as.data.frame <- function(allGs) as.data.frame(matrix(allGs, ncol = 2))

#' @export
fit <- function(x, ...) UseMethod('fit')
#' @export
fit.HQCDP <- function(x, allPars = NULL, method = 'Nelder-Mead') {
	flog.debug('Using method %s', method)
  # here we declare that we want the full output of the rss function for instance
  attr(x, 'complete') <- TRUE
  # get all the parameters to fit if required
  if(is.null(allPars))
    allPars <- getKernelPars(x)
  # reset the bestEvalEnv
  resetBestEval(allPars)
  # the function to internally called by optim
  fn <- function(pars) {
    names(pars) <- names(allPars)
    rssArgs     <- as.list(c(list(x), pars))
    completeVal <- do.call(rss, rssArgs)
    val         <- completeVal$val
    valWeighted <- completeVal$valWeighted
    DoF <- getDoF(x)
    chi2 <- val / DoF
		flog.debug('                        chi2 = %s', round(chi2, 3))
    # store the partial results in the best eval tracker
    saveStep(chi2, val, pars)
		# optimize the log better
    log(valWeighted)
  }

  i <- 1
  while(TRUE) {
    tic()
    bestEval <- get('bestEval', envir = bestEvalEnv)
    lastBestChi2 <- bestEval$chi2
    startPars <- bestEval$pars
    op <- optim(startPars, fn = fn, hessian = FALSE, method = method)
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

getDoF <- function(x, gDeg = 2) {
  # get the experimental points per process
  expPoints <- sum(unlist(lapply(x$processes, function(p) length(expVal(p)))))
  # get the amount of fitting parameters
  # gDeg * number of Reggeons since g(t) = g0 + g1 * t +...
  fitParams <- gDeg * sum(unlist(lapply(x$kernels, function(k) k$numReg)))
  fitParams <- fitParams + length(getKernelPars(x))
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
getSpectra.HQCDP <- function(x, ...) {
	f <- function(t) {
		# we need to initialize the computation on each node
		init()
		val <- list(t = t, spectra = lapply(x$kernels, function(k)
                                  do.call(k$findKernel, list(.t =t, ...))))
		# close the redis connection opened while calling init()
		rredis::redisClose()
		val
	}
  # now compute every kernel for each value of t and the parameters passed
  # this is the most expensive part of the computation and its parallelizable
  # so we use parLapply. Notice that we are already using mclapply in each of the
  # kernel computations so this is the most we can go while nesting parallel computations
  if(!is.null(x$cluster))
    parLapply(x$cluster, getNeededTVals(x), f)
  else
    mclapply(getNeededTVals(x), f, mc.cores = cores)
}
