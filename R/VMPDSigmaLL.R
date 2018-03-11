VMPDSigmaLL <- function(vmName) {
  obs <- DSigma('VMP', 'LL')
  attr(obs, 'vmName') <- vmName
  obs
}

getExternalStateFactor.VMPDSigmaLL <- function(vmpll, Q2 = Q2, alpha = 0) {
  vmName <- attr(vmpll, 'vmName')
  # get the mass of the correspondent vector meson
  mass <- vmpll$data$m[vmpll$data$vmName == vmName][1]
  f1   <- getU1NNMode(Q2 = Q2, alpha = alpha)$dfQ(z)
  f3   <- getU1NormalizableMode(mass = mass, alpha = alpha)$dfQ(z)
  splinefun(z, f1 * f3 / (sqrt(Q2) * mass))
}
