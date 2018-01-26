#' @export
VMPSigma <- function() ProcessObservable('VMPSigma')

#' @export
expVal.VMPSigma <- function(vmp) vmp$data$sigma
#' @export
expErr.VMPSigma <- function(vmp) vmp$data$deltaSigma
#' @export
expKinematics.VMPSigma <- function(vmp) vmp$data[c('m', 'Q2', 'W')]
