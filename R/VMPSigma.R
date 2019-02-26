#' @export
VMPSigma <- function(vmName)
{
  obj <- Sigma('VMP', tmin = -1.2, tmax = 0.0)
  obj$data <- obj$data[obj$data$vmName == vmName,]
  obj
}

#' @export
expVal.VMPSigma <- function(vmp) vmp$data$sigma
#' @export
expErr.VMPSigma <- function(vmp) vmp$data$deltaSigma
#' @export
expKinematics.VMPSigma <- function(vmp) vmp$data[c('m', 'Q2', 'W')]
