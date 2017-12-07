# the data accepted here is the one coming from the kernel function findKernel
#' @export
spectrumBCDQ <- function(data) {
  colors <- c('red', 'blue', 'green', 'cyan')
  trajName <- c('Hard Pomeron', 'Soft Pomeron', '3rd trajectory', '4th trajectory')
  plotSpectrum <- function(indices = NULL, grid = TRUE) {
    if(is.null(indices))
      indices <- 1:length(data)

    z <- get('z', envir = ihqcdEnv)
    plot.new() # reset the plot
    if(grid)
      old.par <- par(mfrow=c(2, 2))
    i <- 0
    plotSet <- FALSE
    info <- c()
    lapply(data, function(d) {
      i <<- i + 1
      if(!is.element(i, indices))
        return()
      # currently in the data comes more than only the spectral data so we have to make sure
      # we plot only the desired part
      if(is.null(d$js))
        return()

      col <- colors[i]
      info <<- c(info, eval(parse(text = paste0('"', trajName[i], ' "~',"j[", i - 1, "]", '~','"="', '~',
                                                format(d$js, digits = if(i > 2) 2 else 3)))))
      if(grid) {
        plot(z, d$u2j, type = 'l', lwd = 3, ylim = c(-20, 15), xlim = c(0, 8),
             ylab = 'U', xlab = 'z', main = info[i])
        abline(h = 0, v = 0, col = "gray60")
        abline(h = c(-5, 5, 10) , v = c(2, 4, 6), col = "lightgray", lty = 3)
        # we need to show explicitly that the potential goes up near the boundary
        if(i > 2) {
          abline(v = 0, lwd = 3)
        }
        # the factor (-1)^(i + 1) is just to match the sign with what is obtained from the data analysis
        lines(d$wf$x, 10 * (1)^(i + 2) * d$wf$y, lwd = 2, lty = 5)
      } else {
        if(!plotSet) {  # first time we have to create a plot, with axis, etc
          plotSet <<- TRUE
          plot(z, d$u2j, type = 'l', lwd = 3, ylim = c(-5, 10), xlim = c(0, 6), col = col, ylab = 'U', xlab = 'z')
          abline(h = c(-5, 5, 10) , v = c(2, 4, 6), col = "lightgray", lty = 3)
        } else  # later we just have to draw lines
          lines(z, d$u2j, lwd = 3, col = col)

        # and always the wavefunctions

        # the factor (-1)^(i + 1) is just to match the sign with what is obtained from the data analysis
        lines(d$wf$x, 5 * (1)^(i + 2) * d$wf$y, lwd = 2, col = col, lty = 5)
      }
    })

    if(!grid) {
      abline(h = 0, v = 0, col = "gray60")
      legend(2, 20, info, col = colors[1: i - 1], lty = c(1, 1), lwd = c(3, 3))
    }
    par(old.par)
  }

  # plot the wave functions of the kernel together with
  # the function obtained acting with the Dperp and Dparallel
  # operators which appear in the non minimal coupling case
  plotNMC <- function(indices = c(1)) {
    As      <- get('As', envir = ihqcdEnv)
    Asder1  <- get('Asder1', envir = ihqcdEnv)
    Asder2  <- get('Asder2', envir = ihqcdEnv)
    Phider1 <- get('Phider1', envir = ihqcdEnv)
    Phider2 <- get('Phider2', envir = ihqcdEnv)
    z       <- get('z', envir = ihqcdEnv)
    # get the computed spectral data
    # and then the specified indices
    actualData <- Filter(function(d) !is.null(d$js), data)[indices]
    # compute Dperp psi
    DperpPsi <- lapply(actualData, function(d) {
      wffun <- splinefun(d$wf$x, d$wf$y)
      wfder1 <- wffun(z, deriv = 1)
      (Asder1 * Phider1 - 0.5 * Asder1^2) * wffun(z) + Asder1 * wfder1
    })
    # compute Dparallel psi
    DparallelPsi <- lapply(actualData, function(d) {
      J <- d$js
      wffun <- splinefun(d$wf$x, d$wf$y)
      wfder1 <- wffun(z, deriv = 1)
      wfder2 <- wffun(z, deriv = 2)
      wfder2 + 2 * (Phider1 - Asder1) * wfder1 +
         (Phider2 + (J - 1.5) * Asder2 - Asder1 * (2 * Phider1 + J - 1) + Phider1^2 + 0.75 * Asder1^2) * wffun(z)
    })

    mapply(function(d, Dp, Dpa, col, i) {
      col <- 'black'
      # the wave function
      plot(d$wf$x, 10 * d$wf$y, type = 'l', lwd = 3, ylim = c(-15, 15), xlim = c(0, 6),
           col = col, ylab = eval(parse(text = paste0('expression(psi[', i - 1, '])'))), xlab = 'z')
      # some grid lines
      abline(h = c(-10, 0, 10) , v = c(2, 4, 6), col = "lightgray", lty = 3)
      # Dperp acting on it
      lines(z, 10 * Dp, lwd = 2, col = col, lty = 3)
      # and Dparallel
      lines(z, 10 * Dpa, lwd = 2, col = col, lty = 2)
    }, actualData, DperpPsi, DparallelPsi, colors[indices], indices)
  }

  s <- list(getData = function() data,
            plot = plotSpectrum,
            plotNMC = plotNMC
            )
  class(s) <- append(class(s), 'spectrumBCDQ')
  s
}

#' @export
plot.spectrumBCDQ <- function(s) s$plot()
