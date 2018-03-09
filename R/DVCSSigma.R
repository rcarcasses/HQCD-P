#' @export
DVCSSigma <- function() Sigma('DVCS', tmin = -1.0, tmax = 0.0)

getExternalStateFactor.DVCSSigma <- getExternalStateFactor.DVCSDSigma

#' @export
enlargeData.DVCSSigma <- function(dvcss, ...) {
  # interesting values of x for producing a plot, for instace.
  Wvals <- seq(20, -2, len = 50)
  Q2s <- sort(unique(dvcss$data$Q2))
  # for each one of the 34 previous values combine them with the xvals
  Reduce(rbind, lapply(Q2s, function(Q2) data.frame(Q2 = rep(Q2, length(Wvals)), W = Wvals)))
}

#' @export
plot.DVCSSigma <- function(dvcss, predicted = NULL) {
  # get all different Q2s
  allQ2s <- sort(unique(dvcss$data$Q2))
  data <- dvcss$data
  values <- data$dsigma
  # prepare the plot
  plot(1e-10, 1,
       log = 'y',
       las = 1,
       xlim = c(30, 165),
       ylim = c(0.07, 40),
       xlab = 'W', ylab = expression(sigma),
       xaxt = 'n', yaxt = 'n')
  hTicks <- 10^(-1:1)
  vTicks <- c(60, 100, 140)
  abline(h = hTicks, v = vTicks, col = 'gray90', lty = 3)
  axis(2, at = hTicks, las = 1, labels = c(expression(10^-1), expression(1), expression(10)))
  axis(1, at = vTicks)
  numCols <- 4
  i <- 1
  cols <- sapply(seq(0.2, 0.7, len = numCols), gray)
  lapply(allQ2s, function(Q2) {
    # get the subset of data for the given Q2
    dataForQ2 <- data[data$Q2 == Q2,]
    colIndex <- mod(i, numCols) + 1
    with(dataForQ2,{
      # draw the experimental points
      lines(W, sigma, type = 'p', pch = 19, col = cols[colIndex], cex = 0.7)
      # draw the experimental error bars
      arrows(W, sigma - deltaSigma, W, sigma + deltaSigma, length = 0.02, angle = 90, code = 3, col = cols[colIndex])
    })
    i <<- i + 1
  })

  invisible(predicted)
}
