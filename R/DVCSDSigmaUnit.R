#' @export
DVCSDSigmaUnit <- function() {
  data <- loadDVCS_dsigma()
  dvcs <- list(data = data)
  class(dvcs) <- append(class(dvcs), 'DVCSDSigmaUnit')
  dvcs
}

#' @export
expVal.DVCSDSigmaUnit <- function(dvcs) dvcs$data$dsigma
#' @export
expErr.DVCSDSigmaUnit <- function(dvcs) dvcs$data$delta_dsigma
#' @export
expKinematics.DVCSDSigmaUnit <- function(dvcs) dvcs$data[c('Q2', 'W', 't')]
