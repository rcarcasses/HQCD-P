#' @export
VMPDSigma <- function(vmName) {
  obs <- DSigma('VMP')
  attr(obs, 'vmName') <- vmName
  # filter experimental data by Vector Meson type
  obs$data <- obs$data[obs$data$vmName == vmName,]
  obs$tt <- VMPDSigmaTT(vmName)
  obs$ll <- VMPDSigmaLL(vmName)
  obs
}

#' @export
predict.VMPDSigma <- function(vmp, Izs, gs, points, ...) {
  amplitude <- getAmplitude(vmp, Izs$tt, gs, points, ...)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / (3 * 16 * pi^2)) * GEVMinus2ToNB * (1 / W^4) * (2*abs(amplitude$tt)^2 + abs(amplitude$ll)^2) / 3
}

#' @export
getCfact.VMPDSigma <- function(vmpds, gs) {
                        CfactName <- paste0('C', attr(vmpds, 'vmName'))
                        if(CfactName %in% names(attributes(gs)))
                          attr(gs, CfactName)
                        else
                          1
                      }
