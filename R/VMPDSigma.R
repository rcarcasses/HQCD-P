#' @export
VMPDSigma <- function(vmName) {
  obs <- DSigma('VMPDSigma')
  attr(obs, 'vmName') <- vmName
  # filter experimental data by Vector Meson type
  obs$data <- obs$data[obs$data$vmName == vmName,]
  obs$tt <- VMPDSigmaTT(vmName)
  obs$ll <- VMPDSigmaLL(vmName)
  obs
}

predict.VMPDSigma <- function(vmp, fns, gs, points, ...) {
  amplitudett <- getAmplitude(vmp$tt, fns, gs, points, ...)
  amplitudell <- getAmplitude(vmp$ll, fns, gs, points, ...)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  GEVMinus2ToNB * (1 / W^4) * (2*abs(amplitudett)^2 + abs(amplitudell)^2) / 3
}
