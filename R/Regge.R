#' @export
regge <- function(model = 'bcde', numReg = 3) {
  p <- potential(model)

  pot <- p$u()
  z   <- p$z
  r   <- NULL

  # set the potential to an external one
  setNewPotential <- function(f, pars) {
    p$setNewPotential(f, pars)
    # update the local potential function
    pot <<- p$u()
  }

  getTrajectories <- function(n = 4, showProgress = TRUE, ...){
    # sometimes extra parameters are passed which are not parameters
    # of the potential function, just ignore them
    fArgs <- getPotentialArgs(pot, ...)
    flog.info('[Regge] getting trajectory for %s', dumpList(fArgs))
    z  <- get('z', envir = ihqcdEnv)
    js <- seq(-0.3, 6.5, len = 100)
    if(showProgress)
      pb <- txtProgressBar(min = min(js), max = max(js), initial = min(js), style = 3)
    t <- function(J) {
      args <- as.list(c(J = J, unlist(fArgs)))
      # cat('calling u2 with args', unlist(args))
      u2j <- do.call(pot, args)
      u2j[1] <- abs(u2j[1]) # make sure the potential goes to plus infinity near z = 0
      dE  <- abs(min(u2j)) / 200;
      # make it at least one for values that give a minimum near 0
      dE <- max(dE, 0.1)
      flush.console()
      setPotential(z, u2j)
      computeSpectrum(n, dE)
      # update progress bar
      if(showProgress)
        setTxtProgressBar(pb, J)
      # return the desired object
      getEnergiesAndIndices()
    }

    getts <- function(J) {
      ts <- t(J)
      ts$energy
    }

    r <<- list(j = js, t = lapply(js, getts))
    r
  }

  plotTrajectories <- function(forceCal = TRUE, n = 4, notJustLines = TRUE, ...) {
    if(is.null(r) || forceCal)
      getTrajectories(n, ...);

    limits <- c(-5, 30)

    for (i in 1:length(r$t[[1]])){
      tr <- list(t = c(), j = c(), j0 = -10, splineInv = c())
      for (k in 1:length(r$j)) {
        # clean all the NAs
        if(is.na(r$t[[k]][i])) next

        tr$t <- c(tr$t, r$t[[k]][i])
        tr$j <- c(tr$j, r$j[[k]])
      }
      # to find the intercepts we find the zeros of the inverse function
      sf <- splinefun(tr$j, tr$t)

      tr$splineInv <- sf
      roots <- uniroot.all(sf, c(min(tr$j), max(tr$j)))
      tr$j0 <- min(roots)
      flog.info(paste('[Regge] j', i - 1,' = ', tr$j0, sep = ''))
      if(i == 1 && notJustLines) {
        plot(tr$t, tr$j, ylim = c(min(r$j), max(r$j)), ylab = 'j(t)', xlab = 't', type = 'n', lwd = 2, xlim = limits)
        abline(h = 0, v = 0, col = "gray10")
        abline(v = 10 * (-1:10), h = c(1:6, seq(0.8, 1.5, by = 0.1)), col = "lightgray", lty = 3)
      }

      lines(tr$t, tr$j, type = 'l', cex = 0.7, xpd = FALSE, lwd = 2)

      labelIndex <- which.min(abs(tr$t - 0.9 * limits[2])) # get the index in the t array with a closest value to the 90% of the right t limit
      boxed.labels(tr$t[labelIndex] , tr$j[labelIndex], labels = format(tr$j0, digits = 3), col = 'blue', cex = 0.7, border = FALSE, bg = 'white')
    }
    pars <- list(...)
    #legend(x = limits[1], y = 4, legend = mapply(function(n, v) getParExpression(n, v), names(pars), pars))
  }

  plotJustLines <- function(n = 4, ...) {
    # call the plot trajectories code asking juts to draw lines
    do.call(plotTrajectories, as.list(c(forceCal = TRUE, n = n, notJustLines = FALSE, list(...))))
  }

  r <- list(get = getTrajectories,
            plot = plotTrajectories,
            plotLines = plotJustLines,
            setNewPotential = setNewPotential,
            ld1 = ld1,
            ld2 = ld2,
            ra = ra,
            wf = wf,
            potential = p)
  class(r) <- append(class(r), 'Regge')
  r
}

#' @export
setNewPotential.Regge <- function(r, ...) r$setNewPotential(...)

# meson masses in GeV

ra <- data.frame(mass = 1e-3 * c(776, 1318, 1689, 1996, 2330, 2450),
                 spin = 1:6)
wf <- data.frame(mass = 1e-3 * c(783, 1275, 1667, 2018, 2250, 2469),
                 spin = 1:6)

# glueball masses in GeV
ld1 <- data.frame(mass  = c(1.475, 2.150, 3.385, 3.640, 4.360), # masses are in GeV, this is the ++ sector
                  spin  = c(0,2,3,4,6),
                  sErr  = c(30, 30, 90, 90, 260),
                  stErr = c(65, 100, 150, 160, 200))
ld1$err <- 1e-3*(ld1$sErr + ld1$sErr)
ld1$err2 <- 2 * ld1$mass * ld1$err

# in the paper the 3++* has a value sqrt(9.44 / 0.9) = 3.239 GeV
# while in the thesis this value doesn't appear and instead a 4++* value apears
# I'll use the thesis's one
ld2 <- data.frame(mass  = c(2.775, 2.880, 0.44 * 10.48), # masses are in GeV, this is the ++ sector
                  spin  = c(0, 2, 4),
                  sErr  = c(70, 100, 0.44 * 5 * 38),
                  stErr = c(120, 130, 130))
ld2$err <- 1e-3*(ld2$sErr + ld2$sErr)
ld2$err2 <- 2 * ld2$mass * ld2$err

# by now we don't want the scalar glueball data, neither spin 3, which should belong to the odderon trajectory
ld1 <- ld1[-c(1,3),]
ld2 <- ld2[-c(1),]

#' @export
plot.reggeBCDQ <- function(r, ...) r$plot(forceCal = TRUE, n = 4,...)

#' @export
drawMesonsMasses <- function() {
  lines(ra$mass^2, ra$spin, type = 'p', pch = 5, col = 'green')
  lines(wf$mass^2, wf$spin, type = 'p', pch = 6, col = 'green')
}

#' @export
drawGlueballMasses <- function() {
  lines(ld1$mass^2, ld1$spin, type = 'p', pch = 20)
  # plot the error bars
  arrows(ld1$mass^2 - ld1$err2, ld1$spin,
         ld1$mass^2 + ld1$err2, ld1$spin,
         length = 0.05, angle = 90, code = 3)

  lines(ld2$mass^2, ld2$spin, type = 'p', pch = 21)
  # plot the error bars
  arrows(ld2$mass^2 - ld2$err2, ld2$spin,
         ld2$mass^2 + ld2$err2, ld2$spin,
         length = 0.05, angle = 90, code = 3)

}
