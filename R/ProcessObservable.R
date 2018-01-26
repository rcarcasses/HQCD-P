#' @export
ProcessObservable <- function(className) {
  obs <- list()
  class(obs) <- append(class(obs), className)
  data <- loadData(obs)
  obs <- c(obs, list(data = data))
  class(obs) <- append(class(obs), className)
  obs
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
