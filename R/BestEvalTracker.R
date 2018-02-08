
#' This environment contain the best values for the fit
#' in real time
#' @export
bestEvalEnv <- new.env()
assign('startGs', NULL, envir = bestEvalEnv)

#' @export
saveStep <- function(chi2, val, pars) {
  evalStep <- list(chi2 = chi2, val = val, pars = pars)
  assign('lastEval', evalStep, envir = bestEvalEnv)
  # check if this is the best evaluation so far
  bestEval <- get('bestEval', envir = bestEvalEnv)
  if(bestEval$val > val) {
    assign('bestEval', evalStep, envir = bestEvalEnv)
    str <- paste('-> chi2 =', round(chi2, 4),' pars =', do.call(paste, as.list(round(pars, 3))))
    flog.debug(str)
  }
}

#' @export
resetBestEval <- function(pars) assign('bestEval', list(val = 1e50, pars = pars), envir = bestEvalEnv)
