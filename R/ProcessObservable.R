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
predict.ProcessObservable <- function(obs, spectra, points = NULL, fns = NULL, gs = NULL, ...) {
  # if no points are provided, use the experimental ones
  if(is.null(points))
    points <- expKinematics(obs)
  if(is.null(fns))
    fns <- getFns(obs, points, spectra)
  if(is.null(gs))
    gs <- getBestGs(obs, fns)
  # call the next function
  arg <- list(...)
  arg$points  <- points     # add the experimental points of this process
  arg$fns     <- fns        # add the fns of this process
  arg$gs      <- gs         # add the experimental points of this process
  arg$generic <- 'predict'  # set the next method to call
  arg$object  <- obs        # and pass the same object as argument
  do.call(NextMethod, arg)  # finally call the next method in the chain
}

#' Set the non minimal coupling attribute which is used to known
#' if the non minimal coupling contribution will be considered or not
#' @return same object with the NMC attribute set
#' @export
setNMC <- function(x, value = TRUE) {
  attr(x, 'NMC') <- value
  x
}
#' Indicate whether this process will be predicted considering a non
#' minimal coupling or not
#' @return TRUE or FALSE
#' @export
useNMC <- function(x) attr(x, 'NMC')

#' Every experiment have a set of values for t for which
#' the kernels have to be computed. This information should
#' be returned by this function
#' @export
getNeededTVals <- function(x) UseMethod('getNeededTVals')
#' @export
getNeededTVals.default <- function(x) 'getNeededTVals should be custom implemented for each ProcessObservable subtype'
#' @export
getFns <- function(x, ...) UseMethod('getFns')
#' @export
getFns.default <- function(x) paste('getFns has to be implemented for this object with classes', class(x))

#' Given an ProcessObservable object and its associated fns already computed
#' find the best coefficients, gs, such that the sum of difference squared
#' weigthed by the inverse of the error squared (rss) between the prediction
#' and the experimental values of the process is minimized.
#' @export
getBestGs <- function(x, ...) UseMethod('getBestGs')
#' @export
getBestGs.default <- function(x) paste('getBestGs has to be implemented for this object with classes', class(x))

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
fNNMC <- function(x, ...) UseMethod('fNNMC')

#' @export
getJs <- function(x, ...) UseMethod('getJs')
getJs.default <- function(x) 'getJs called in an object with no implementation'

#' @export
plotSpectrum <- function(x, ...) UseMethod('plotSpectrum')
plotSpectrum.default <- function(x) 'plotSpectrum called in an object with no implementation'
