#' @export
DVCSDSigma <- function() DSigma('DVCSDSigma')

getExternalStateFactor.DVCSDSigma <- function(dvcs, Q2 = Q2, ...) {
  getU1NNMode(Q2 = Q2)$fQ
}

