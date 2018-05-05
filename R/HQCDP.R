#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function(alpha = 0, hOrder = 2, rootRejectionWeight = 1, rootRejectionCutoff = 0.02) {
  h <- list(processes = list(), kernels = list())
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
  # add the constraint for the intercept of the soft pomeron
  # the value of this attribute will be used as weight while fitting
  attr(h, 'addSPconstraint')     <- 1e6
  attr(h, 'hOrder')              <- hOrder
  attr(h, 'alpha')               <- alpha
  attr(h, 'rootRejectionWeight') <- rootRejectionWeight
  attr(h, 'rootRejectionCutoff') <- rootRejectionCutoff
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
addProcessObservable.HQCDP <- function(p, ...) {
  Reduce(function(h, x) {
    # insert the non minimal coupling attribute in the child processes
    attr(x, 'alpha') <- attr(h, 'alpha')
    if(is.element('ProcessObservable', class(x)))
      h$processes <- append(h$processes, list(x))
    else
      flog.error(paste('Adding an object of unknown classes:', class(x)))
    # return the modified HQCD object
    h
  }, list(...), init = p)
}

#' @export
print.HQCDP <- function(x) {
  cat('HQCDP:\n * processes:', unlist(lapply(x$processes, function(proc) tail(class(proc), n = 1))), '\n * kernels:', length(x$kernels), '\n')
}

#' This computes the sum of the rss for all the processes configured.
#' This is the function to minimize in the optimization procedure.
#' @export
rss.HQCDP <- function(x, pars, zstar, hpars) {
  # warn the user that there are configurations
  # that are still required
  if(length(x$processes) == 0 || length(x$kernels) == 0) {
    cat('Please complete your HQCDP configuration first\n')
    return()
  }
  # first we need to compute the spectra of the kernels
	# this is a parallelized call
  spectra <- getSpectra(x, pars)
  # now we find all the Izs for each one of the processes
  allProcIzs    <- lapply(x$processes, getIzs, spectra = spectra)
  allProcIzsBar <- lapply(x$processes, getIzsBar, spectra = spectra, zstar = zstar, hpars = hpars)

  # find the rss for each process for the given values of gs and fns
  val <- sum(unlist(mapply(function(proc, Izs, IzsBar) {
          rss(proc, Izs = Izs, IzsBar = IzsBar)
        }, x$processes, allProcIzs, allProcIzsBar)))
  # add SP constraint here
  valWeighted <- 0
  if(!is.null(attr(x, 'addSPconstraint'))) {
    spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
    jsSP <- spectraForTZero[[1]][[2]]$js
    # cat('SP intercept', jsSP, '\n')
    valWeighted <- val + attr(x, 'addSPconstraint') * (jsSP - 1.09)^2
  }
  # find all roots and push them away from the zeros
  js <- getJs(x, spectra)
  roots <- uniroot.all(Vectorize(function(J) H(J, hpars)), c(0.8, 1.2))
  valWeighted <- valWeighted + sum(
    unlist(
      lapply(js, function(J) {
        # check if the intercept is close enough to a root of H(J)
        # and add a penalization if so
        attr(x, 'rootRejectionWeight') * sum(unlist(lapply(roots, function(root) {
          if(abs(J - root) < attr(x, 'rootRejectionCutoff')) 1 else 0
        })))
  })))

  # return the result in the appropiated format
  if(is.null(attr(x, 'complete')))
    val
  else
    list(val = val, valWeighted = valWeighted)
}

#' @export
fit <- function(x, ...) UseMethod('fit')
#' @export
fit.HQCDP <- function(x, pars = NULL, zstar = 0.6, hpars = NULL, method = 'Nelder-Mead') {
	flog.debug('Using method %s, number of cores: %s', method, cores)
  DoF <- getDoF(x)
	flog.debug('number of degrees of freedom: %s', DoF)
  # here we declare that we want the full output of the rss function for instance
  attr(x, 'complete') <- TRUE
  # get all the parameters to fit if required
  if(is.null(pars))
    pars <- getKernelPars(x)
  if(is.null(hpars))
    hpars <- c(-166.51, 171.41, -201.52)
  # reset the bestEvalEnv
  initPars <- c(pars, zstar, hpars)
  cat('init par', initPars, '\n')
  resetEvalTracker(initPars)
  # the function to internally called by optim
  fn <- function(fitPars) {
    mark <- length(fitPars) - (attr(x, 'hOrder') + 2)
    if(mark != 0) {
      rssParsIndices <- 1:mark
      rssPars <- fitPars[rssParsIndices]
      zstar   <- fitPars[-rssParsIndices][1]
      hpars   <- fitPars[-rssParsIndices][-1]  # remove the rss pars
      names(rssPars) <- names(pars)        # put the right names
    } else {
      zstar   <- fitPars[1]
      hpars   <- fitPars[-1]
      rssPars <- c()
    }
    #cat('rssPars', rssPars, ' zstar', zstar, 'hpars', hpars, '\n')
    lastEval    <- get('lastEval', envir = bestEvalEnv)
    completeVal <- tryStack(
      rss(x, pars = rssPars, zstar = zstar, hpars = hpars)
    )
    # if the computation ends with an error a string is returned with its description
    if(is.character(completeVal))
      return(1e3 * (1 + 0.05 * runif(1)))
    # otherwise complete the computation and save it
		val         <- completeVal$val
    valWeighted <- completeVal$valWeighted
    chi2        <- val / DoF
    flog.debug('          pars = %s', do.call(paste, as.list(round(fitPars, 6))))
		flog.debug('          chi2 = %s', round(chi2, 3))
    # store the partial results in the best eval tracker
    saveStep(chi2, val, fitPars)
		# optimize the log better
    log(valWeighted)
  }

  tic()
  bestEval <- get('bestEval', envir = bestEvalEnv)
  lastBestChi2 <- bestEval$chi2
  startPars    <- bestEval$pars
  op <- optim(startPars,
							fn = fn,
							hessian = FALSE, method = method)
  bestEval <- get('bestEval', envir = bestEvalEnv)
  newBestChi2 <- bestEval$chi2
  if(abs(newBestChi2 - lastBestChi2) < 1e-3) # the new iteration wasn't better than the one before
    break

  lastBestChi2 <- newBestChi2
  exectime <- toc()
  flog.debug(paste('[HQCDP]  -::- fit completed in', seconds_to_period(round(exectime$toc - exectime$tic)), ' min, best chi2', bestEval$chi2))
  op
}

getJs.HQCDP <- function(x, spectra) {
  # get the t = 0 spectra
  unlist(lapply(Filter(function(s) s$t == 0, spectra)[[1]]$spectra, function(spec) lapply(spec, `[[`, 'js')))
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
  fitParams <- attr(x, 'hOrder') + 2 + length(getKernelPars(x))
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
plot.HQCDP <- function(x, predicted = NULL, pars = NULL, gs = NULL, dry = FALSE) {
  if(is.null(predicted)) {
    # we need to compute the spectra for an enlarged set of t values
    # convert gs to a data frame, if required
    if(!is.data.frame(gs))
      gs <- gs.as.data.frame(x, gs)
    # get the plot points
    plotPoints <- enlargeData(x)
    # compute the spectrum, now with the particular t values needed
    ts <- unique(unlist(lapply(p$processes, function(proc) {
      if('Sigma' %in% class(proc))
        enlargeKinematicsWithTs(proc)$t
      else
        enlargeData(proc)$t
    })))
    pb <- txtProgressBar(min = 0, max = 100, initial = 1, style = 3)
    spectra <- getSpectra(x, pars, ts)
    cat('\n') # put the progress bar in a new line
    setTxtProgressBar(pb, 30)
    # get the Izs for the plot points
    i <- 0
    allProcIzs <- mapply(function(proc, points) {
      setTxtProgressBar(pb, 30 + 30 * i / length(x$processes))
      i <<- i + 1
      list(getIzs(proc, spectra = spectra, points = points))
    }, x$processes, plotPoints)
    # find the predictions for the plot points for each one of the processes
    i <- 0
    predicted <- mapply(function(proc, procIzs, procPlotPoints) {
      pred <- predict(proc, points = procPlotPoints, Izs = procIzs, gs = gs)
      setTxtProgressBar(pb, 60 + 40 * i / length(x$processes))
      i <<- i + 1
      #cat('predictions found', unlist(pred), '\n')
      list(cbind(procPlotPoints, data.frame(predicted = pred)))
    }, x$processes, allProcIzs, plotPoints)
  }

  # call the plot function on each on
  if(!dry)
    mapply(plot, x$processes, predicted)
  # return invisibly the computed points, useful to speed up future plots
  invisible(predicted)
}
