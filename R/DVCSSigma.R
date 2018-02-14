#' @export
DVCSSigma <- function() Sigma('DVCSSigma')

getExternalU1wf.DVCSSigma <- function(dvcs, Q2 = Q2, ...) {
  getU1NNMode(Q2 = Q2)$fQ
}
