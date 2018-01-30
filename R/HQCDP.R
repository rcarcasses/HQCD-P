#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function() {
  h <- list(processes = list(), kernels = list())
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
  # add the constraint for the intercept of the soft pomeron
  # the value of this attribute will be used as weight while fitting
  attr(h, 'addSPconstraint') <- 1e6
  h
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
  # first we need to compute the spectra of the kernels
  spectra <- getSpectra(x, ...)
  val <- sum(unlist(lapply(x$processes, rss, spectra = spectra)))
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
    list(val = val, valWeighted = valWeighted)
}

#' @export
fit <- function(x, ...) UseMethod('fit')
#' @export
fit.HQCDP <- function(x, allPars = NULL) {
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
    # store the partial results in the best eval tracker
    saveStep(chi2, val, pars)
    valWeighted
  }

  i <- 1
  while(TRUE) {
    tic()
    bestEval <- get('bestEval', envir = bestEvalEnv)
    lastBestChi2 <- bestEval$chi2
    startPars <- bestEval$pars
    op <- optim(startPars, fn = fn, hessian = FALSE)
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

getDoF <- function(x){
  # get the experimental points per process
  expPoints <- sum(unlist(lapply(x$processes, function(p) length(expVal(p)))))
  # get the amount of fitting parameters
  fitParams <- sum(unlist(lapply(x$kernels, function(k) k$numReg)))
  fitParams <- fitParams + length(getKernelPars(x))
  expPoints - fitParams
}
#' @export
getSpectra <- function(x, ...) UseMethod('getSpectra')
#' @export
getSpectra.default <- function(x) 'getSpectra not defined in current object'
#' @export
getSpectra.HQCDP <- function(x, ...) {
  # now compute every kernel for each value of t and the parameters passed
  lapply(getNeededTVals(x), function(t)
                              list(t = t,
                                   spectra = lapply(x$kernels, function(k)
                                     do.call(k$findKernel, list(.t =t, ...)))))
}
