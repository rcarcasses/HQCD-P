#' @export
VMPDSigma <- function() ProcessObservable('VMPDSigma')

#' @export
expVal.VMPDSigma <- function(vmp) vmp$data$dsigma
#' @export
expErr.VMPDSigma <- function(vmp) vmp$data$deltaDSigma
#' @export
expKinematics.VMPDSigma <- function(vmp) vmp$data[c('m', 'Q2', 'W', 't')]
