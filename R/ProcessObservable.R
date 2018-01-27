#' @export
ProcessObservable <- function(className) {
  obs <- list()
  class(obs) <- append(class(obs), className)
  data <- loadData(obs)
  obs <- c(obs, list(data = data))
  class(obs) <- append(class(obs), c('ProcessObservable', className))
  obs
}

#' @export
predict.ProcessObservable <- function(obs, points = NULL, ...) {
  # if no points are provided, use the experimental ones
  if(is.null(points))
    points <- expKinematics(obs)
  # call the next function
  NextMethod('predict', obs, points, ...)
}

#' @export
rssComplete <- function(x, ...) UseMethod('rssComplete')
#' @export
rssComplete.ProcessObservable <- function(obs, ...) {
  pred <- predict(obs, ...)
  exp  <- expVal(obs)
  err  <- expErr(obs)
  list(rss = sum(((pred - exp) / err)^2))
}

#' @export
expVal <- function(x) UseMethod('expVal')
#' @export
expErr <- function(x) UseMethod('expErr')
#' @export
expKinematics <- function(x) UseMethod('expKinematics')
#' @export
loadData <- function(x) UseMethod('loadData')

#' @export
getAllFns <- function(x, ...) UseMethod('getAllFns')
getAllFns.default <- function(x, ...) flog.error('getAllFns called in an object with no implementation')

#' @export
getFns <- function(x, ...) UseMethod('getFns')
getFns.default <- function(x) 'getFns called in an object with no implementation'

#' @export
setNewPotential <- function(x, ...) UseMethod('setNewPotential')
setNewPotential.default <- function(x) 'setNewPotential called in an object with no implementation'

#' @export
getJs <- function(x, ...) UseMethod('getJs')
getJs.default <- function(x) 'getJs called in an object with no implementation'

#' @export
plotSprectrum <- function(x, ...) UseMethod('plotSprectrum')
plotSprectrum.default <- function(x) 'plotSprectrum called in an object with no implementation'
