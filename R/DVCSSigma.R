#' @export
DVCSSigma <- function() ProcessObservable('DVCSSigma')

#' @export
expVal.DVCSSigma <- function(dvcs) dvcs$data$sigma
#' @export
expErr.DVCSSigma <- function(dvcs) dvcs$data$deltaSigma
#' @export
expKinematics.DVCSSigma <- function(dvcs) dvcs$data[c('Q2', 'W')]

