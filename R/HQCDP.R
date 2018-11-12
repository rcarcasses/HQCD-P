#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function(alpha = 0,
                  fixed = list(),
                  rsslog = FALSE,
                  rootRejectionWeight = 1, rootRejectionCutoff = 0.02,
                  H = NULL, hparsInitDefault = NULL) {
  h <- list(processes = list(), kernels = list())
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
  if(is.null(hparsInitDefault))
    hparsInitDefault <- c(-1.677950, 1.714511, -2.015598, 0, 0)
  # if H was not passed, use a default one
  if(is.null(H)) {
    flog.warn('Using default H(J) function (are you sure?)')
    H <- function(J, hpars)
           exp(100 * (hpars[1] + hpars[2] * J + hpars[3] * log(J) + hpars[4] * J * log(J) + hpars[5] * J^2))
  }

  # add the constraint for the intercept of the soft pomeron
  # the value of this attribute will be used as weight while fitting
  attr(h, 'addSPconstraint')     <- 1e4
  attr(h, 'fixed')               <- fixed
  attr(h, 'alpha')               <- alpha
  attr(h, 'rootRejectionWeight') <- rootRejectionWeight
  attr(h, 'rootRejectionCutoff') <- rootRejectionCutoff
  attr(h, 'hparsInitDefault')    <- hparsInitDefault
  attr(h, 'cacheSpectra')        <- FALSE
  attr(h, 'H')                   <- H
  attr(h, 'rsslog')              <- rsslog
  attr(h, 'showSpectraProgress') <- TRUE
  attr(h, 'saveLastSpectra')     <- TRUE
  attr(h, 'useTVals')            <- c()
  # compute the gns
  h
}


#' Given an actual list or number, replace the NA
#' on it with the values in thegiven 'field of the
#'  attribute 'fixed' of the x HQCDP object passed
#' @export
completeWithFixedVal <- function(x, actual, field) {
  # if there are not any NA, then there is nothing to be replaced
  if(!anyNA(actual))
    return(actual)
  # get the replacement(s) for the NA, if any
  replacement <- attr(x, 'fixed')[[field]]
  if(is.null(replacement)) {
    flog.warn('There is no replacements in field %s of attribute fixed of the HQCDP object for the NA values', field)
    return(actual)
  }
  # store where there are not NA in replacement, to be replaced
  naPlaces <- !is.na(replacement)
  # if is a vector replace the NA with the correspondent
  # values from the replacement vector, else just replace the number
  if(is.vector(actual))
    actual[naPlaces] <- replacement[naPlaces]
  else
    actual <- replacement

  if(anyNA(actual))
    flog.warn('The replacement using the field %s of the attribute fixed of the HQCDP object seems to be incomplete', field)

  actual
}

#' @export
addFixedValues <- function(x, actual, NAIndices, field) {
  if(length(NAIndices) == 0)
    return(actual)
  # add the missing NA
  for (i in 1:length(NAIndices)) {
    actual <- append(actual, NA, after = NAIndices[i] - 1)
  }
  # return the competed value
  completeWithFixedVal(x, actual, field)
}


#' @export
getNeededTVals.HQCDP <- function(x) {
    # find the sorted unique union of all the processes needed t values
    sort(unique(c(0, unlist(lapply(x$processes, getNeededTVals)))))
}
#' @export
enlargeData.HQCDP <- function(x) lapply(x$processes, enlargeData)

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
    # copy attributes that need to propagate
    attr(x, 'alpha')  <- attr(h, 'alpha')
    attr(x, 'rsslog') <- attr(h, 'rsslog')
    attr(x, 'H')      <- attr(h, 'H')
    # take care if there is a child dsigma object
    if(!is.null(x$dsigma)) {
      attr(x$dsigma, 'alpha')  <- attr(h, 'alpha')
      attr(x$dsigma, 'rsslog') <- attr(h, 'rsslog')
      attr(x$dsigma, 'H')      <- attr(h, 'H')
    }

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
  # find all roots and push them away from the poles
  js <- getJs(x, spectra)
  roots <- uniroot.all(Vectorize(function(J) attr(x, 'H')(J, hpars)), c(0.8, 1.2))
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
fit.HQCDP <- function(x, pars = NULL, zstar = 0.565, hpars = NULL, method = 'Nelder-Mead') {
	flog.debug('Using method %s, number of cores: %s', method, cores)
  # here we declare that we want the full output of the rss function for instance
  attr(x, 'complete') <- TRUE
  # get all the parameters to fit if required
  if(is.null(pars))
    pars <- getKernelPars(x)

  if(is.null(hpars))
    hpars <- attr(x, 'hparsInitDefault')

  # we need to save the positions where the NA are, which
  # means such parameters are fixed with the values used in the
  # constructor of x
  parsNAIndices  <- which(is.na(pars))
  isZstarNA      <- is.na(zstar)
  hparsNAIndices <- which(is.na(hpars))

  # remove all the NAs
  initPars <- c(
          if(length(parsNAIndices) > 0)  pars[-parsNAIndices]   else pars,
          if(isZstarNA)                  NULL                   else zstar,
          if(length(hparsNAIndices) > 0) hpars[-hparsNAIndices] else hpars
        )
  # reset the bestEvalEnv
  resetEvalTracker(initPars)
  DoF <- getDoF(x, initPars)
  flog.debug('init pars %s', do.call(paste, as.list(round(initPars, 3))))
	flog.debug('number of degrees of freedom: %s', DoF)
  # the function to internally called by optim
  fn <- function(fitPars) {
    mark <- length(getKernelPars(x)) - length(parsNAIndices)
    if(mark != 0) {
      # complete the full parameters demanded by the rss function if NA appears
      rssParsIndices <- 1:mark
      rssPars <- addFixedValues(x, fitPars[rssParsIndices], parsNAIndices, 'pars')
      zstar   <- if(isZstarNA) completeWithFixedVal(x, NA, 'zstar') else fitPars[-rssParsIndices][1]
      hpars   <- if(isZstarNA) addFixedValues(x, fitPars[-rssParsIndices], hparsNAIndices, 'hpars') else addFixedValues(x, fitPars[-rssParsIndices][-1], hparsNAIndices, 'hpars')  # remove the rss pars
      names(rssPars) <- names(pars)            # put the right names
    } else {
      zstar   <- if(isZstarNA) completeWithFixedVal(x, NA, 'zstar') else fitPars[1]
      hpars   <- if(isZstarNA) addFixedValues(x, fitPars, hparsNAIndices, 'hpars') else addFixedValues(x, fitPars[-1], hparsNAIndices, 'hpars')
      rssPars <- c()
    }
    #cat('rssPars', rssPars, ' zstar', zstar, 'hpars', hpars, '\n')
    lastEval    <- get('lastEval', envir = bestEvalEnv)
    completeVal <- #tryStack(
      rss(x, pars = rssPars, zstar = zstar, hpars = hpars)
    #)
    # if the computation ends with an error a string is returned with its description
    if(is.character(completeVal) || is.na(completeVal$val))
      return(1e3 * (1 + 0.05 * runif(1)))
    # otherwise complete the computation and save it
		val         <- completeVal$val
    valWeighted <- completeVal$valWeighted
    chi2        <- val / DoF
    flog.debug('          pars = %s', do.call(paste, as.list(round(fitPars, 6))))
		flog.debug('          chi2 = %s', round(chi2, 3))
    # store the partial results in the best eval tracker
    saveStep(chi2, val, fitPars)
		# optimize the log better ?
    valWeighted
  }

  tic()
  bestEval <- get('bestEval', envir = bestEvalEnv)
  lastBestChi2 <- bestEval$chi2
  startPars    <- bestEval$pars
  op <- optim(startPars,
							fn = fn,
							hessian = FALSE, method = method, control = list(maxit = 10000))
  bestEval <- get('bestEval', envir = bestEvalEnv)
  newBestChi2 <- bestEval$chi2
  if(abs(newBestChi2 - lastBestChi2) < 1e-3) # the new iteration wasn't better than the one before
    break

  lastBestChi2 <- newBestChi2
  exectime <- toc()
  flog.debug(paste('[HQCDP]  -::- fit completed in', seconds_to_period(round(exectime$toc - exectime$tic)), ' min, best chi2', bestEval$chi2))
  op
}

#' @export
getJs.HQCDP <- function(x, spectra) {
  # get the t = 0 spectra
  unlist(lapply(Filter(function(s) s$t == 0, spectra)[[1]]$spectra, function(spec) lapply(spec, `[[`, 'js')))
}

getKernelPars <- function(x) {
  allPars <- unlist(lapply(x$kernels, function(k) k$optimPars))
  keep    <- !duplicated(names(allPars))
  allPars[keep]
}

#' @export
getDoF <- function(x, pars) {
  # get the experimental points per process
  expPoints <- sum(unlist(lapply(x$processes, function(p) length(expVal(p)))))
  expPoints -  length(pars)
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
#' Get the spectra of a configured HQCDP object.
#' @param x a HQCDP object, already configured.
#' @param pars additional parameters which will be passed
#' down while computing the spectrum.
#' @param ts the values of t for which the spectrum will be computed. If this parameter is not passed
#' then the t values needed are deduced from the union of the t values needed in the computation of each
#' one of the process of the HQCDP object. If such object, contains an attribute 'useTVals' being a non
#' zero length vector, then the values in this vector are used to compute a spectra which later is used
#' to find the spectra for the real t values needed, as defined by the union of the t values needed of
#' all the processes, by using interpolation.
#' @export
getSpectra.HQCDP <- function(x, pars = NULL, ts = NULL) {
  # if the attribute 'useTVals' have been set with a length bigger than zero,
  # then use those, in this case we have to later interpolate the values obtained
  if(length(attr(x, 'useTVals')) > 0)
    ts <- attr(x, 'useTVals')
  else if(is.null(ts))
    ts <- getNeededTVals(x)

  # in order to speed up computations we cache, if enabled the option,
  # the computation of the spectra as well
  # create a unique key for this function call arguments values
  key <- digest(c(pars, ts), algo = 'md5')
  if(attr(x, 'cacheSpectra')) {
    cached <- rredis::redisGet(key)
    if(!is.null(cached)) {
      #flog.trace('Using cached spectra list')
      return(cached)
    } else
      flog.trace('Spectra will be computed and cached')
  }

	numRegs <- unlist(lapply(x$kernels, function(k) k$numReg))
	# unwrap the computation of each one of the Reggeon data for each one of the kernels
	unwrappedFunCalls <- Reduce(function(acct, t) {
  	as.list(c(acct,
  	  Reduce(function(acc, k) {
    	  as.list(c(acc,
    	    lapply(1:k$numReg, function(n) {
    	      # it may be the case that for a given kernel we are passing
    	      # more parameters than those needed by it, the findReggeonData
    	      # function will appropriately discard them
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
  spectra <- convertRawSpectra(
    if(Sys.info()['sysname'] == 'Linux') {
      # it is a good idea apparently to explicitly split the function calls
      # into chunks, otherwise it seems like a large unwrappedFunCalls list
      # causes memory issues
      # see the chunks example at https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
      chunks <- split(unwrappedFunCalls, cut(seq_along(unwrappedFunCalls),
                             ceiling(length(unwrappedFunCalls) / (2 * cores)), labels = FALSE))
      if(attr(x, 'showSpectraProgress'))
        pb <- progress_bar$new(format = " computing spectra [:bar] :percent eta: :eta",
                                total = length(chunks), clear = FALSE, width= 60)
      # flatten list one level only
      unlist(lapply(chunks, function(chunk) {
        # update the progress bar
        if(attr(x, 'showSpectraProgress'))
          pb$tick()
        # do the batch computation to find kernel pieces
        mclapply(chunk, function(f) f(), mc.cores = cores)
      }), recursive = FALSE)
      #mclapply(unwrappedFunCalls, function(f) f(), mc.cores = cores)
    }
    else {
      if(attr(x, 'showSpectraProgress'))
        pb <- progress_bar$new(format = " computing spectra [:bar] :percent eta: :eta",
                               total = length(unwrappedFunCalls), clear = FALSE, width= 60)
      lapply(unwrappedFunCalls, function(f) {
        # update the progress bar
        if(attr(x, 'showSpectraProgress'))
          pb$tick()
        # call the function and return its value
        f()
      })
    }
  , numRegs, ts)

  # if the attribute 'useTVals' have been set with a length bigger than zero,
  # then here we need to use the computed spectra for the 'useTVals' values
  # and interpolate it for the real t values needed
  if(length(attr(x, 'useTVals')) > 0) {
    # get the real ones needed
    ts <- getNeededTVals(x)
    flog.trace('Interpolating spectra from %s points to %s', length(attr(x, 'useTVals')), length(ts))
    spectra <- interpolateSpectraTo(spectra, ts)
  }

  # cache if enabled and needed
  if(attr(x, 'cacheSpectra'))
    if(is.null(rredis::redisGet(key))) {
      flog.trace('Saving full spectra to cache')
      rredis::redisSet(key, spectra)
    }
  # store the last spectra for debug purposes
  if(attr(x, 'saveLastSpectra'))
    saveRDS(spectra, file = '~/lastSpectra.rds')
  # return the computed object
  spectra
}

interpolateSpectraTo <- function(spectra, ts) {
  # get all the t for which the spectra was computed
  usedTs <- unlist(lapply(spectra, `[[`, 't'))
  # get the x values for the wave functions used
  xs <- spectra[[1]]$spectra[[1]][[1]]$wf$x
  # flatten different kernels reggeons and put them in the same level
  flattenSpectra <- unlist(lapply(spectra, `[[`, 'spectra'), recursive = FALSE)
  # get the set of indices
  indices <- 1:length(flattenSpectra[[1]])

  # the index here is the index of the reggeon when one takes into account all the
  # kernels, like for example with a two kernel configuration of two reggeons
  # index 3 is the first reggeon of the second kernel
  interpolateReggeon <- function(index) {
    # get the values of J for each t and its derivative
    Jt <- unlist(lapply(lapply(flattenSpectra, `[[`, index), `[[`, 'js'))
    dJdt <- unlist(lapply(lapply(flattenSpectra, `[[`, index), `[[`, 'dJdt'))
    # create the interpolation function for the J values dependence with t
    Jfun <- splinefun(usedTs, Jt)
    dJdtfun <- splinefun(usedTs, dJdt)
    # get the values of wf$y for each t
    ys <- lapply(lapply(lapply(flattenSpectra, `[[`, index), `[[`, 'wf'), `[[`, 'y')
    # for each position x produce an interpolation function that describe
    # how the point i of the wave functio index changes with t
    ysfun <- lapply(1:length(xs), function(i) {
      # get all the values of y for the given i index, which correspond to xs[i]
      yi <- unlist(lapply(ys, `[[`, i))
      # make a spline and return it
      splinefun(usedTs, yi)
    })
    list(ysfun = ysfun, Jfun = Jfun, dJdtfun = dJdtfun, index = index)
  }
  # do the interpolation per index
  spectraInterpolationPerIndex <- lapply(indices, interpolateReggeon)
  # now with the previous computation done we need to actually use the interpolation
  # functions found to find the approximated spectra for the actually needed
  # values of t
  pb <- progress_bar$new(format = " interpolating spectra [:bar] :percent eta: :eta",
                          total = length(ts), clear = FALSE, width= 60)
  lapply(ts, function(tval) {
    computedSpectra <- lapply(indices, function(index) {
      sindex <- spectraInterpolationPerIndex[[index]]
      # get the structure of the given index object, we pick the first element
      # of flattenSpectra since the structure is the same for any of it
      obj <- flattenSpectra[[1]][[index]]
      # update the relevant information
      obj$js   <- sindex$Jfun(tval)
      obj$wf$y <- unlist(lapply(sindex$ysfun, function(f) f(tval)))
      obj$dJdt <- sindex$dJdtfun(tval)
      # finally return the modified object, notice that properties
      # as name, index and numReg are not touched
      obj
    })
    # now the previous computation has the structure of a flattened spectra, we need
    # to put it with the right structure taking into account the amount of
    # actual kernels
    computedSpectraStructured <- list(list())
    kerIndex <- 1
    regIndex <- 1
    for (i in indices) {
      computedSpectraStructured[[kerIndex]][[regIndex]] <- computedSpectra[[i]]
      # if this is the data for the last reggeon of the present kernel
      # then update the kernel index, recall that we are assumming here that
      # the spectra list is organized such that each reggeon index property
      # matches its position inside the list
      if(i != indices[length(indices)] && computedSpectra[[i]]$index == computedSpectra[[i]]$numReg) {
        kerIndex <- kerIndex + 1
        computedSpectraStructured[[kerIndex]] <- list()
        regIndex <- 0
      }
      regIndex <- regIndex + 1
    }

    # finaly return the t and the correspondent spectra computed
    list(t = tval, spectra = computedSpectraStructured)
  })
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
plot.HQCDP <- function(x, predicted = NULL, pars = NULL, zstar = NULL, hpars = NULL, dry = FALSE) {
  if(is.null(predicted)) {
    # we need to compute the spectra for an enlarged set of t values
    # get the plot points
    plotPoints <- enlargeData(x)
    # compute the spectrum, now with the particular t values needed
    ts <- unique(unlist(lapply(p$processes, function(proc) {
      if('Sigma' %in% class(proc))
        enlargeKinematicsWithTs(proc)$t
      else
        enlargeData(proc)$t
    })))
    flog.trace('Number of t values %s', length(ts))
    spectra <- getSpectra(x, pars = pars, ts = ts)
    pb <- progress_bar$new(format = " computing Izs [:bar] :percent eta: :eta",
                            total = length(x$processes), clear = FALSE, width= 60)
    # get the Izs for the plot points
    allProcIzs <- mapply(function(proc, points) {
      pb$tick()
      list(getIzs(proc, spectra = spectra, points = points))
    }, x$processes, plotPoints)
    pb <- progress_bar$new(format = " computing IzsBar [:bar] :percent eta: :eta",
                           total = length(x$processes), clear = FALSE, width= 60)
    # get the IzBars for the plot points
    allProcIzsBar <- mapply(function(proc, points) {
      pb$tick()
      list(getIzsBar(proc, spectra = spectra, points = points, zstar = zstar, hpars = hpars))
    }, x$processes, plotPoints)

    pb <- progress_bar$new(format = " predicting values for plot points [:bar] :percent eta: :eta",
                           total = length(x$processes), clear = FALSE, width= 60)
    # find the predictions for the plot points for each one of the processes
    predicted <- mapply(function(proc, Izs, IzsBar, points) {
      pred <- predict(proc, Izs = Izs, IzsBar = IzsBar, points = points)
      pb$tick()
      #cat('predictions found', unlist(pred), '\n')
      list(cbind(points, data.frame(predicted = pred)))
    }, x$processes, allProcIzs, allProcIzsBar, plotPoints)
  }

  # call the plot function on each on
  if(!dry)
    mapply(plot, x$processes, predicted)
  # return invisibly the computed points, useful to speed up future plots
  invisible(predicted)
}
