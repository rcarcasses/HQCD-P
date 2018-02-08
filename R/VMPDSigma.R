#' @export
VMPDSigma <- function(vmName) {
  obs <- DSigma('VMPDSigma')
  attr(obs, 'vmName') <- vmName
  # filter experimental data by Vector Meson type
  obs$data <- subset(obs$data, vmName = vmName)
  obs
}

getExternalU1wf.DVCSDSigma <- function(vmp, Q2 = Q2, ...) {
  vmName <- attr(vmp, 'vmName')
  #TODO: find the right expression
  getU1NNMode(Q2 = Q2)$fQ
}

