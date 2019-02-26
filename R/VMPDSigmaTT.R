VMPDSigmaTT <- function(vmName) {
  obs <- DSigma('VMP', 'TT')
  attr(obs, 'vmName') <- vmName
  obs
}

getExternalStateFactor.VMPDSigmaTT <- function(vmptt, Q2 = Q2, alpha = 0) {
  f1   <- getU1NNMode(Q2 = Q2, alpha = alpha)$fQ(z)
  f3   <- getU1NormalizableMode()$fm(z)
  splinefun(z, f1 * f3)
}

getCfact.VMPDSigmaTT <- getCfact.VMPDSigma
