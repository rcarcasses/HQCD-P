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
addKernel <- function(x, ...) UseMethod('addKernel')
#' @export
addKernel.default <- function(x, ...) 'calling addKernel in the wrong object'
#' @export
addKernel.HQCDP <- function(h, x) {
  if(is.element('kernelUnit', class(x)))
    h$kernels <- append(h$kernels, list(x))
  else
    flog.error(paste('Adding an object of unknown classes:', class(x)))
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
getSpectra <- function(x, ...) UseMethod('getSpectra')
#' @export
getSpectra.default <- function(x) 'getSpectra not defined in current object'
#' @export
getSpectra.HQCDP <- function(x, ...) {
  # first get a list of values of t for which the kernels are needed
  tVals <- unlist(lapply(x$processes, function(p) getNeededTVals(p)))
  # now compute every kernel for each value of t and the parameters passed
  lapply(tVals, function(t) list(t = t, kernels = lapply(x$kernels, function(k) do.call(k$findKernel, list(...)))))

}
