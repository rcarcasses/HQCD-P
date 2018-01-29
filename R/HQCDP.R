#' This is the main function of this library.
#' It allow to define a model with many kernels that can be
#' tested again a configurable list of experimental obsevables
#' @export
HQCDP <- function() {
  h <- list(processes = list(), kernels = list())
  class(h) <- c('HQCDP', class(h))  # pay attention to the class name
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
  sum(unlist(lapply(x$processes, rss, spectra = spectra)))
}

#' @export
fit <- function(x, ...) UseMethod('fit')
#' @export
fit.HQCDP <- function(x) {
  # get all the parameters to fit
  allPars <- getKernelPars(x)
  fn <- function(pars) {
    names(pars) <- names(allPars)
    rssArgs <- as.list(c(list(x), pars))
    val <- do.call(rss, rssArgs)
    DoF <- getDoF(x)
    chi2 <- val / DoF
    cat( ' chi2 =', chi2, '\n')
    val
  }
  op <- optim(allPars, fn = fn, hessian = FALSE)
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
