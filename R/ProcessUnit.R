#' @export
getAllFns <- function(x) UseMethod('getAllFns')
getAllFns.default <- function(x) flog.error('getAllFns called in an object with no implementation')

#' @export
getFns <- function(x) UseMethod('getFns')
getFns.default <- function(x) flog.error('getFns called in an object with no implementation')

#' @export
setNewPotential <- function(x) UseMethod('setNewPotential')
setNewPotential.default <- function(x) flog.error('setNewPotential called in an object with no implementation')

#' @export
getJs <- function(x) UseMethod('getJs')
getJs.default <- function(x) flog.error('getJs called in an object with no implementation')

#' @export
plotSprectrum <- function(x) UseMethod('plotSprectrum')
plotSprectrum.default <- function(x) flog.error('plotSprectrum called in an object with no implementation')

#' @export
plotGlueballMasses <- function(x) UseMethod('plotGlueballMasses')
plotGlueballMasses.default <- function(x) flog.error('plotGlueballMasses called in an object with no implementation')

#' @export
plotMesonsMasses <- function(x) UseMethod('plotMesonsMasses')
plotMesonsMasses.default <- function(x) flog.error('plotMesonsMasses called in an object with no implementation')
