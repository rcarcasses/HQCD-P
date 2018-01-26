#' @export
DVCSSigmaUnit <- function() {
  data <- loadData.DVCSSigma()
  dvcs <- list(data = data)
  class(dvcs) <- append(class(dvcs), 'DVCSSigmaUnit')
  dvcs
}

#' @export
expVal.DVCSSigmaUnit <- function(dvcs) dvcs$data$sigma
#' @export
expErr.DVCSSigmaUnit <- function(dvcs) dvcs$data$deltaSigma
#' @export
expKinematics.DVCSSigmaUnit <- function(dvcs) dvcs$data[c('Q2', 'W')]
