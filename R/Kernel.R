#' this is a stateful object: its content can be configured per run, this is
#' the amount of kernels and its related amount of reggeons, etc
#' @export
Kernel <- function() {
  kernels <- list()

  addKernel <- function(ker) {
    kernels <- c(kernels, ker)
  }

  getKernels <- function() kernels

  k <- list(addKernel = addKernel,
            getKernels = getKernels)
  k <- class(Kernel) <- append(class(k), 'KernelHQCDP')  # pay attention to the class name
  k
}

#' @export
addKernel <- function(x) UseMethod('addKernel')
#' @export
addKernel.default <- function(x) 'addKernel not defined in current object'
#' @export
addKernel.KernelHQCDP <- function(ks, k) ks$addKernel(k)

#' @export
getKernels <- function(x) UseMethod('getKernels')
#' @export
getKernels.default <- function(x) 'getKernels not defined in current object'
#' @export
getKernels.KernelHQCDP <- function(ks) ks$getKernels(k)

#' @export
getSpectra <- function(x, ...) UseMethod('getSpectra')
#' @export
getSpectra.default <- function(x) 'getSpectra not defined in current object'
#' @export
getSpectra.KernelHQCDP <- function(ks, ...) {
  lapply(getKernels(ks), function(k) do.call(k$findKernel, list(...)))
}
