#' @export
DVCSSigmaUnit <- function() {
  data <- loadDVCS_sigma()
  dvcs <- list(data = data)
  class(dvcs) <- append(class(dvcs), 'DVCSSigmaUnit')
  dvcs
}

#' @export
expVal.DVCSSigmaUnit <- function(dvcs) dvcs$data$sigma
#' @export
expErr.DVCSSigmaUnit <- function(dvcs) dvcs$data$delta_sigma
#' @export
expKinematics.DVCSSigmaUnit <- function(dvcs) dvcs$data[c('Q2', 'W')]
