#' @export
Sigma <- function(procName, tmin, tmax, ...) {
  # add the generic DSigma class
  obs <- ProcessObservable(paste0(procName, 'Sigma'))
  class(obs) <- append(class(obs), 'Sigma', after = 2)
  attr(obs, 'tmin') <- tmin
  attr(obs, 'tmax') <- tmax
  # create an internal dsigma object and
  # enhance artificially the data including N values of t for
  # each entry point. These will be required for the computation
  # of the integral
  obs$dsigma <- DSigma(procName, ...)
  obs
}

#' For these processes this enlargement process will insert
#' ficticious values of t for each one of the kinematical entries
#' such that the differential cross-section can be computed
#' @export
enlargeKinematicsWithTs <- function(sigma, data = NULL) {
  if(is.null(data))
    data <- expKinematics(sigma)
  ts <- getNeededTVals(sigma)
  as.data.frame(lapply(as.data.frame(Reduce(function(acc, r) {
    val <- rbind(acc, r)
    rownames(val) <- NULL
    val
  }, unlist(
    apply(data, 1,
          function(row)
            lapply(ts, function(t) c(list(t = t), row))),
    recursive = FALSE)
  )), `mode<-`, 'numeric'))
}
# Here we need to take the t values from t = -1 to t = 0
#' @export
getNeededTVals.Sigma <- function(x) seq(attr(x, 'tmin'), attr(x, 'tmax'), len = 20)

#' Predicts the values of F2 for the points passed
#' @param sigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param gs A list of list which contains the coefficients of  g(t) = g0 + g1 * t + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.Sigma <- function(sig, Izs, IzsBar, points, ...) {
  # compute all the needed differential cross-sections
  # For each value of Q2 and W we need to insert many different values of t
  sig$dsigma$data <- enlargeKinematicsWithTs(sig, points)
  dsigma <- predict(sig$dsigma, Izs = Izs, IzsBar = IzsBar, ...)
  ts <- getNeededTVals(sig)
  numBlocks <- length(points[[1]])
  #flog.trace('[Sigma] number of t blocks %s', numBlocks)
  # we have created above blocks of row in the dataframe that correspond
  # to the same values of Q2 and W with different values of t between -1 and 0
  # now we need to compute the differential cross-section for each value of t
  # inside a given block, and integrate it to get the total cross-section
  # associated with the given block, or equivalently, to the values of Q2 and W
  sigma <- unlist(lapply(1:numBlocks, function(i) {
    # get the dsigma data computed for the value of Q2 and W, etc...
    # in the present block
    ds <- dsigma[1:length(ts) + (i - 1) * length(ts)]
    # check if there is any NA and remove it
    na.pos <- which(is.na(ds))
    if(length(na.pos) > 0) {
      flog.warn('NA values found at Sigma computation, discarding them...')
      ds <- ds[-na.pos]
      ts <- ts[-na.pos]
    }
    # if there are not enough points we will return a big number and throw an error
    if(length(ds) < 2) {
      flog.error('Not enough information to predict Sigma')
      return(1e20)
    }
    # we now need to compute the cross section as the integral over [-1,0] of dsigma
    dsFun <- splinefun(ts, ds)
    #flog.trace('[Sigma] computing integral from %s to %s', attr(sig, 'tmin'), attr(sig, 'tmax'))
    integral <- integrate(dsFun, attr(sig, 'tmin'), attr(sig, 'tmax'), stop.on.error = FALSE)
    integral$value
  }))
  sigma
}

#' @export
getIzs.Sigma <- function(sigma, points, spectra) {
  # if a custom way of computing the integral was defined, then notify the child dsigma object about it
  if('IzN' %in% names(attributes(sigma)))
    attr(sigma$dsigma, 'IzN') <- attr(sigma, 'IzN')
  # return the z integrals for the enhanced points of the correspondent dsigma object
  getIzs(sigma$dsigma, spectra = spectra, points = enlargeKinematicsWithTs(sigma, points))
}

#' @export
getIzsBar.Sigma <- function(sigma, points, spectra, zstar, hpars) {
  # if a custom way of computing the integral was defined, then notify the child dsigma object about it
  if('IzNBar' %in% names(attributes(sigma)))
    attr(sigma$dsigma, 'IzNBar') <- attr(sigma, 'IzNBar')
  # return the z bar integrals for the enhanced points of the correspondent dsigma object
  getIzsBar(sigma$dsigma, spectra = spectra, points = enlargeKinematicsWithTs(sigma, points), zstar = zstar, hpars = hpars)
}

#' @export
expVal.Sigma <- function(s) s$data$sigma
#' @export
expErr.Sigma <- function(s) s$data$deltaSigma
#' @export
expKinematics.Sigma <- function(s) s$data[c('Q2', 'W')]

