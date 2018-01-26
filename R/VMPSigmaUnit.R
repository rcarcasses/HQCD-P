#' @export
VMPSigmaUnit <- function() {
  data <- loadVMP_dsigma()
  vmp <- list(data = data)
  class(vmp) <- append(class(vmp), 'VMPSigmaUnit')
  vmp
}

#' @export
expVal.VMPSigmaUnit <- function(vmp) vmp$data$sigma
#' @export
expErr.VMPSigmaUnit <- function(vmp) vmp$data$delta_sigma
#' @export
expKinematics.VMPSigmaUnit <- function(vmp) vmp$data[c('m', 'Q2', 'W')]
