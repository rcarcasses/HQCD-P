#' @export
DVCSDSigma <- function() DSigma('DVCS')

getExternalStateFactor.DVCSDSigma <- function(dvcs, Q2 = Q2, ...) {
  getU1NNMode(Q2 = Q2)$fQ
}

#' predict the values
#' @param dsigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param hs A data frame which contains the coefficients of  h(J) = h0 + h1 * (J-1) + ...
#' for each one of the Reggeons.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.DVCSDSigma <- function(dsigma, Izs, IzsBar, points, alpha = 0, ...) {
  amplitude <- getAmplitude(dsigma, Izs, IzsBar, points, ...)
  absAmplitudeSquared <- if (abs(alpha) >  0) {
    amplitudeNMC1 <- getAmplitudeNMC1(dsigma, Izs, IzsBar, points, ...)
    amplitudeNMC2 <- getAmplitudeNMC2(dsigma, Izs, IzsBar, points, ...)
    # TODO: return the differential cross sections
    ReAMC <- Re(amplitude)
    ImAMC <- Im(amplitude)
    ReANMC1 <- Re(amplitudeNMC1)
    ImANMC1 <- Im(amplitudeNMC1)
    ReANMC2 <- Re(amplitudeNMC2)
    ImANMC2 <- Im(amplitudeNMC2)
    ImAMC^2 + 2 * ImAMC * ImANMC1 + ImANMC1^2 + points$t * ImAMC * ImANMC2 +
      points$t * ImANMC1 * ImANMC2 + (3/8) * points$t^2 * ImANMC2^2 + ReAMC^2 +
      2 * ReAMC * ReANMC1 + ReANMC1^2 + points$t * ReAMC * ReANMC2 +
      points$t * ReANMC1 * ReANMC2 + (3/8) * points$t^2 * ReANMC2^2
  } else {
    abs(amplitude)^2
  }
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / (16 * pi^2)) * GEVMinus2ToNB * (1 / W^4) * absAmplitudeSquared
}

#' @export
plot.DVCSDSigma <- function(dvcsds, predicted, numGraphs = 4) {
  plot.new()
  par(mfrow = c(2, 2), las = 1)
  # get all different Q2s
  allQ2s <- sort(unique(dvcsds$data$Q2))
  # do as many plot as demanded
  invisible(lapply(1:numGraphs, function(n) {
    Q2s <- allQ2s[sapply(0:(length(allQ2s) - 1), mod, numGraphs) + 1 == n]
    data <- dvcsds$data[dvcsds$data$Q2 %in% Q2s,]
    values <- data$dsigma
    # prepare the plot
    plot(0,
         log = 'y',
         las = 1,
         xlim = c(0, 1),
         ylim = c(0.02, 50),
         xlab = '-t', ylab = expression(d * sigma / dt),
         xaxt = 'n', yaxt = 'n')
    hTicks <- seq(0, 1, len = 6)
    vTicks <- c(0.1, 1, 10)
    abline(v = hTicks, h = vTicks, col = 'gray90', lty = 3)
    axis(1, at = hTicks)
    axis(2, at = vTicks)
    i <- 1
    cols <-  sapply(seq(0, 1, len = 5), gray)
    legendItems <- list()
    lapply(Q2s, function(Q2) {
      # get the data for this Q2
      dataQ2 <- data[data$Q2 == Q2,]
      # get all the different W for this Q2
      Ws <- unique(dataQ2$W)
      cat('Ws', Ws, '\n')
      j <- 1
      lapply(Ws, function(W) {
        dataQ2andW <- dataQ2[dataQ2$W == W,]
        size <-  0.5 * j + 0.5
        pch <- 15 + i
        col <- cols[j]
        with(dataQ2andW, {
          # plot the data
          lines(-t, dsigma, type = 'p', pch = pch, col = col)
          # plot the errors
          arrows(-t, dsigma - deltaDSigma, -t, dsigma + deltaDSigma, length = 0.02 * j, angle = 90, code = 3, col = col)
          # find the predicted values and plot them
          predictedQ2andW <- predicted[predicted$Q2 == Q2[1] & predicted$W == W[1],]
          # order the entries
          predictedQ2andW <- predictedQ2andW[order(predictedQ2andW$t),]
          lines(-predictedQ2andW$t, predictedQ2andW$predicted, col = col)
        })
        legendItems[[length(legendItems) + 1]] <<- list(pch = pch, col = col, W = W, Q2 = Q2)
        # get the predicted values for the given values of W and Q2
        j <<- j + 1
      })
      i <<- i + 1
    })
    # add the legend
    txt <- sapply(legendItems, function(item) as.expression(substitute(Q^2 == Q2 ~~ 'W' == W, item)))
    legend(0.6, 70, txt, pch = sapply(legendItems, `[[`, 'pch'), col = sapply(legendItems, `[[`, 'col'), bty = 'o', y.intersp = 0.9, cex = 0.8, box.lwd = 0, box.col = "white",bg = "white")
  }))
  invisible(par())
}
