#' @export
VMPDSigma <- function(vmName) {
  obs <- DSigma('VMPDSigma')
  attr(obs, 'vmName') <- vmName
  # filter experimental data by Vector Meson type
  obs$data <- subset(obs$data, vmName = vmName)
  obs
}

predict.VMPDSigma <- function(vmp, fns, gs, points, ...) {
  # compute g(t) for each corresponding fn, this return a dataframe
  # were each column is the value (a vector) of g(t) for the given values
  # of t and the n column is the g(t) of the n fn in the fns dataframe
  # see the test file for some explanation 'tests/testthat/test_DSigma.R'
  gts <- apply(gs, 1, function(row) {
    rowSums(t(row * t(outer(t, 0:(length(gs) - 1), `^`))))
  })
  # compute the amplitude
  amplitude <- rowSums(fns * gts, na.rm = TRUE)
  # get the Ws
  W <- points$W
  # return the differential cross sections
  (1 / W^4) * abs(amplitude)^2
}
