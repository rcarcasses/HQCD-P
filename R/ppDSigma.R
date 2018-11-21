#' @export
ppDSigma <- function() DSigma('pp')

#' @export
getNeededTVals.pp <- function(pp) unique(pp$data$minus_t)

#' predict the values
#' @param dsigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.ppDSigma <- function(dsigma, Izs, IzsBar, points, ...) {
  amplitude <- getAmplitude(dsigma, Izs, IzsBar, points, ...)
  # get the Ws
  W <- points$W
  # return the differential cross sections, converted to MB
  df <- GEVMinus2ToMB * (1 / (16 * pi * W^4)) * abs(amplitude)^2
  #dput(df)
  df
}

plot.ppDSigma <- function(pp, p, pars, hpars, zstar) {
  points <- expKinematics(pp)
  ts <- expKinematics(pp)$t
  plot.new()
  #par(mfrow = c(2, 2), las = 1, mar = c(4, 4, 1, 2))
  data <- pp$data
  values <- data$dsigma
  # prepare the plot
  plot(0,
       log = 'y',
       las = 1,
       xlim = c(0, 4),
       ylim = c(1e-7, 60),
       xlab = '-t', ylab = expression(d * sigma / dt),
       xaxt = 'n', yaxt = 'n')
  hTicks <- seq(0, 4, len = 6)
  vTicks <- seq(1e-6,10, len = 10)
  abline(v = hTicks, h = vTicks, col = 'gray90', lty = 3)
  axis(1, at = hTicks)
  axis(2, at = vTicks)
  # plot the data
  points(-data$t, data$dsigma, type = 'p', pch = 1, col = "black")
  # find the predicted
  spectra <- getSpectra(p, pars = pars, ts = ts)
  pb <- progress_bar$new(format = " computing Izs [:bar] :percent eta: :eta",
                         total = length(p$processes), clear = FALSE, width= 60)
  izs <- getIzs(pp, spectra = spectra)
  izbars <- getIzsBar(pp, spectra = spectra, zstar = zstar, hpars = hpars)
  # compute the prediction
  prediction <- predict(pp, Izs = izs, IzsBar = izbars, points = points)
  points(-data$t, prediction, col = "blue")
}
