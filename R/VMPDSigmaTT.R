#' @export
VMPDSigmaTT <- function(vmName) {
  obs <- DSigma('VMPDSigmaTT')
  attr(obs, 'vmName') <- vmName
  obs
}

getExternalStateFactor.VMPDSigmaTT <- function(vmptt, Q2 = Q2, alpha = 0) {
  vmName <- attr(vmp, 'vmName')
  # get the mass of the correspondent vector meson
  mass <- vmp$data$m[vmp$data$vmName == vmName][1]
  f1   <- getU1NNMode(Q2 = Q2, alpha = alpha)$fQ(z)
  f3   <- getU1NormalizableMode(mass = mass, alpha = alpha)$fQ(z)
  splinefun(z, f1 * f3)
}
