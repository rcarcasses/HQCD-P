#' @export
VMPDSigmaTT <- function(vmName) {
  obs <- DSigma('VMPDSigmaTT')
}

getExternalU1wf.VMPDSigmaTT <- function(vmptt, Q2 = Q2, ...) {
  vmName <- attr(vmp, 'vmName')
  #TODO: find the right expression
  wf <- getU1NNMode(Q2 = Q2)$fQ(z)*
}
