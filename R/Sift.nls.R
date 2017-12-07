#' Perform a nls minimization with a custom function.
#' @param data Data should contain a "y" and a "err" columns, the rest are considered as independent variables.
#' @param options It is possible to change the name of the dependent variable and the error by modifying this argument.
#' @export
sift.nls <- function(fn = NULL, data = NULL, start = NULL, traces = TRUE, optimMethod = 'Nelder-Mead', chi2max = 6, options = list(yName = 'y', errName = 'err', normal = FALSE)) {

  yName   <- options$yName
  errName <- options$errName

  if(is.null(data[[yName]]) || is.null(data[[yName]])) {
    cat('\nPlease make sure you are providing the right names for the dependent variable and error columns\n')
    return()
  }

  dof <- length(data[[yName]]) - length(start)

  chi2vecFun <- function (alpha) {
    by(data, 1:nrow(data), function(r) {
                          # a list with the remaining row values
                          x <- r
                          x[[yName]]   <- NULL
                          x[[errName]] <- NULL

                          chi2(alpha, x, r[[yName]], r[[errName]])
                       })
  }

  # the common version is like this
  rss <- function(...) {
    chi2vec <- chi2vecFun(...)
    dof     <- length(data[[yName]]) - length(...)
    s       <- sum(chi2vec)
    chi2    <- s / dof
    list(rss = s, dof = dof, chi2 = chi2)
  }

  # here we define the "robust" formula to minimize which gives the same result as minimizing the squared deviations
  # later I should workout a vectorized version
  Lambda <-  function(...) {
    chi2vec    <- chi2vecFun(...)
    sLambda <- sum(unlist(lapply(chi2vec, function(c2) log(1 + 0.18 * c2))))
    s       <- sum(chi2vec)
    chi2    <- s / dof
    list(rss = sLambda, dof = dof, chi2 = chi2, rssNormal = s)
  }

  chi2 <- function(alpha, x, y, err) {
    yt <- do.call(fn, c(alpha, x))
    ((yt - y) / err)^2
  }

  optimFun <- Lambda
  if(options$normal)
    optimFun <- rss

  bestEval <- NULL

  optimFunWrapper <- function(...) {
    s   <- optimFun(...)
    c2  <- do.call(chi2tot, list(...))
    if(traces)
      cat('\n', names(start), '\n', ..., ' chi2 =>', c2, '\n')

    # keep track of the best result found
    if(!is.na(s)) {
      if(is.null(bestEval))
        bestEval <<- list(pars = list(...), val = s, chi2 = c2)
      else if(bestEval$val > s )
        bestEval <<- list(pars = list(...), val = s, chi2 = c2)
    }
    s
  }

  #op <- optim(start,
  #            fn      = optimFunWrapper,
  #            method  = optimMethod)

  # just a dumb fix
  # bestEval$pars <- bestEval$pars[[1]]
  sift.data <- function(chi2max, ...) {
    chi2vec   <- chi2vecFun(list(...))
    newData   <- Filter(function(c2) c2 < chi2max, chi2vec)
    sieveData <- Filter(function(c2) c2 > chi2max, chi2vec)
    s         <- sum(chi2vec)
    sNew      <- sum(newData)
    chi2      <- s / dof
    chi2new   <- sNew / dof
    list(rss = s, dof = dof, chi2 = chi2, chi2new = chi2new, data = newData, removed = sieveData)
  }

  #return(bestEval)
  s <- list( optimFun = optimFun,
            sift.data = sift.data
            )
  class(s) <- append(class(s), 'sift.nls')
  s
}
