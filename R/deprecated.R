
# deprecated
#' @export
computeWFs <- function(b0 = 260, b1 = 260, c0 = 0, c1 = 0, ls = 0.112) {
  u2j1 <- U2Jextra(j1, b0, b1, c0, c1, ls)
  z   <- get('z', envir = bcdqEnv)
  dE  <- abs(min(u2j1)) / 200;
  # make it at least one for values that give a minimum near 0
  dE <- max(dE, 1)
  bcdqLog("j1 = ", j1, " dE1 =", dE, "\n")
  setPotential(z[1:length(u2j1)], u2j1)
  computeSpectrum(1, dE)
  e1  <- getEnergies()[1]
  wf1 <- getWavefunctions()[[1]]

  u2j2 <- U2Jextra(j2, b0, b1, c0, c1, ls)
  dE  <- abs(min(u2j2)) / 200;
  dE <- max(dE, 1)
  bcdqLog("j2 = ", j2, " dE2 =", dE, "\n")
  setPotential(z[1:length(u2j2)], u2j2)
  computeSpectrum(2, dE)
  e2  <- getEnergies()[2]
  wf2 <- getWavefunctions()[[2]]
  bcdqLog('e1 =', e1, ' e2 =', e2, '\n')

  assign('u2j1', u2j1, envir = bcdqEnv)
  assign('u2j2', u2j2, envir = bcdqEnv)
  assign('wf1', wf1, envir = bcdqEnv)
  assign('wf2', wf2, envir = bcdqEnv)
}

#' @export
plotKernelWFs <- function(create = TRUE, max1 = 50, max2 = 50, ...) {
  u2j1 <- get('u2j1', envir = bcdqEnv)
  u2j2 <- get('u2j2', envir = bcdqEnv)
  wf1 <- get('wf1', envir = bcdqEnv)
  wf2 <- get('wf2', envir = bcdqEnv)
  z  <- get('z', envir = bcdqEnv)

  wf1$y = (max1 / max(wf1$y)) * wf1$y
  wf2$y = (max2 / max(wf2$y)) * wf2$y

  if(create)
    plot.new()

  options = list(...)
  #do.call(plot, append(list(z, u2j1, type = 'l', lwd = 3, col = 'blue'), options))
  lines(z, u2j1, type = 'l', lwd = 3, col = 'red', ylim = c(-50, 100), xlim = c(0, 2.5))
  lines(z, u2j2, type = 'l', lwd = 3, col = 'blue')
  lines(wf1$x, wf1$y, col = 'red')
  lines(wf2$x, wf2$y, col = 'blue')
}

#' @export
quantumAnalysisBCDQ <- function(b0 = 260, b1 = 260, c0 = 0, c1 = 0, ls = 0.112) {
  j0 <- get("j0", envir = heraEnv)
  u2j <- U2Jextra(j0[1], b0, b1, c0, c1, ls)
  dE  <- abs(min(u2j)) / 200;
  bcdqLog('dE =', dE)
  z  <- get('z', envir = bcdqEnv)
  quantumAnalysis(z[1:length(u2j)], u2j, 2, 400, dE, FALSE, ylim = c(-50, 50), xlim = c(0, 4))
  cat(getEnergies())
}
