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
