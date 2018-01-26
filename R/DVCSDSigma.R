#' @export
DVCSDSigma <- function() ProcessObservable('DVCSDSigma')

#' @export
expVal.DVCSDSigma <- function(dvcs) dvcs$data$dsigma
#' @export
expErr.DVCSDSigma <- function(dvcs) dvcs$data$deltaDSigma
#' @export
expKinematics.DVCSDSigma <- function(dvcs) dvcs$data[c('Q2', 'W', 't')]
