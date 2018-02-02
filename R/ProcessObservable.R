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
  arg <- list(...)
  arg$points  <- points
  arg$generic <- 'predict'
  arg$object  <- obs
  do.call(NextMethod, arg)
}

#' Every experiment have a set of values for t for which
#' the kernels have to be computed. This information should
#' be returned by this function
#' @export
getNeededTVals <- function(x) UseMethod('getNeededTVals')
#' @export
getNeededTVals.default <- function(x) 'getNeededTVals should be custom implemented for each ProcessObservable subtype'

#' @export
rss <- function(x, ...) UseMethod('rss')
#' @export
rss.ProcessObservable <- function(obs, ...) {
  pred <- predict(obs, ...)
  # recall that predict give additional information if the attribute 'complete' is set
  if(!is.null(attr(obs, 'complete')))
    pred <- pred$val
  exp  <- expVal(obs)
  err  <- expErr(obs)
  val  <- sum(((pred - exp) / err)^2)
  # return the full information according with the attribute 'complete' of the object
  if(is.null(attr(obs, 'complete')))
    val
  else
    list(val = val)
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
setNewPotential <- function(x, ...) UseMethod('setNewPotential')
setNewPotential.default <- function(x) 'setNewPotential called in an object with no implementation'

#' @export
fN <- function(x, ...) UseMethod('fN')

#' @export
computeAmplitude <- function(x, ...) UseMethod('computeAmplitude')

#' @export
getJs <- function(x, ...) UseMethod('getJs')
getJs.default <- function(x) 'getJs called in an object with no implementation'

#' @export
plotSprectrum <- function(x, ...) UseMethod('plotSprectrum')
plotSprectrum.default <- function(x) 'plotSprectrum called in an object with no implementation'
