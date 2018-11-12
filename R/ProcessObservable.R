ProcessObservable <- function(className) {
  obs <- list()
  class(obs) <- append(class(obs), className)
  data <- loadData(obs)
  obs <- c(obs, list(data = data))
  class(obs) <- append(class(obs), c('ProcessObservable', className))
  attr(obs, 'alpha')  <- 0
  attr(obs, 'rsslog') <- FALSE
  obs
}

#' @export
predict.ProcessObservable <- function(obs,
                                      spectra = NULL,
                                      points = NULL,
                                      Izs = NULL,
                                      IzsBar = NULL, ...) {
  # if there is not enough data complaint
  if(is.null(spectra) && is.null(Izs)) {
    flog.error(paste('Not enough information to make a prediction over the object with classes', class(obs)))
    return()
  }
  # if no points are provided, use the experimental ones
  if(is.null(points))
    points <- expKinematics(obs)
  if(is.null(Izs))
    Izs <- getIzs(obs, points = points, spectra = spectra)
  if(is.null(IzsBar))
    IzsBar <- getIzsBar(obs, points = points, spectra = spectra)

  # call the next function
  arg <- list(...)
  arg$points  <- points     # add the experimental points of this process
  arg$Izs     <- Izs        # add the Izs of this process
  arg$IzsBar  <- IzsBar     # add the IzsBar of this process
  arg$generic <- 'predict'  # set the next method to call
  arg$object  <- obs        # and pass the same object as argument
  do.call(NextMethod, arg)  # finally call the next method in the chain
}

#' @export
getIzs.ProcessObservable <- function(obs, spectra = NULL, ...) {
  # if there is not enough data complaint
  if(is.null(spectra)) {
    flog.error(paste('Not enough information to find the z integrals over the object with classes', class(obs)))
    return()
  }
  do.call(NextMethod, injectPoints(obs, 'getIzs', spectra = spectra, ...)) # finally call the next method in the chain
}

#' @export
getIzsBar.ProcessObservable <- function(obs, spectra = NULL, ...) {
  # if there is not enough data complaint
  if(is.null(spectra)) {
    flog.error(paste('Not enough information to find the z bar integrals over the object with classes', class(obs)))
    return()
  }
  do.call(NextMethod, injectPoints(obs, 'getIzsBar', spectra = spectra, ...)) # finally call the next method in the chain
}

injectPoints <- function(obs, genName, ...) {
  arg <- list(...)
  # if no points are provided, use the experimental ones
  if(is.null(arg$points))
    arg$points <- expKinematics(obs)

  arg$generic <- genName # set the next method to call
  arg$object  <- obs     # and pass the same object as argument
  arg
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
getIzs <- function(x, ...) UseMethod('getIzs')
#' @export
getIzs.default <- function(x) paste('getIzs has to be implemented for this object with classes', class(x))
#' @export
getIzsBar <- function(x, ...) UseMethod('getIzsBar')
getIzsBar.default <- function(x) paste('getIzsBar has to be implemented for this object with classes', class(x))

#' @export
IzNBar <- function(x, ...) UseMethod('IzNBar')
#' @export
IzNBar.ProcessObservable <- function(obs, kin, spec, zstar, hpars) {
  J     <- spec$js
  wf    <- spec$wf
  dJdt  <- spec$dJdt
  t1fun <- splinefun(z, exp((-J + 0.5) * As + Phi))
  t2fun <- getExternalProtonFactor()
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)$value
  # get the function H from the observable attribute H
  H <- attr(obs, 'H')
  # if it doesn't exist use a default one
  # THIS FUNCTION DEPENDING IF WE CONSIDER PP SCATTERING OR NOT. THINK ABOUT THIS LATER
  if(is.null(H))
    H <- function(J, hpars)
           exp( hpars[1] + hpars[2] * (J - 1) + hpars[3] * (J - 1)^2 )
  # Return IzNBar
  iznbar <- H(J, hpars) * (1i + 1 / tan( 0.5 * pi * J) ) * dJdt * integral
  return(iznbar)
}

#' @export
rss <- function(x, ...) UseMethod('rss')
#' @export
rss.ProcessObservable <- function(obs, ...) {
  val  <- sum((diffObsWeighted(obs, ...))^2)
  # return the full information according with the attribute 'complete' of the object
  if(is.null(attr(obs, 'complete')))
    val
  else
    list(val = val)
  }


diffObsWeighted <- function(x, ...) UseMethod('diffObsWeighted')
diffObsWeighted.ProcessObservable <- function(obs, ...) {
  Opred <- predict(obs, ...)
  # recall that predict give additional information if the attribute 'complete' is set
  if(!is.null(attr(obs, 'complete')))
    Opred <- Opred$val

  Oexp  <- expVal(obs)
  Oerr  <- expErr(obs)

  if(attr(obs, 'rsslog'))
    Oexp * log(Opred / Oexp) / Oerr
  else
    (Opred - Oexp) / Oerr
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
IzN <- function(x, ...) UseMethod('IzN')

#' @export
IzNNMC <- function(x, ...) UseMethod('IzNNMC')

#' @export
getJs <- function(x, ...) UseMethod('getJs')
getJs.default <- function(x) 'getJs called in an object with no implementation'

#' @export
plotSpectrum <- function(x, ...) UseMethod('plotSpectrum')
plotSpectrum.default <- function(x) 'plotSpectrum called in an object with no implementation'

#' @export
enlargeData <- function(x, ...) UseMethod('enlargeData')
enlargeData.default <- function(x, ...) 'enlargeData called in an object with no implementation'
