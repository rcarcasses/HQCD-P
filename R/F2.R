#' @export
F2 <- function() ProcessObservable('F2')

#' For DIS we only need t = 0 value
#' @export
getNeededTVals.F2 <- function(f2) c(0)

#' Predicts the values of F2 for the points passed
#' @param f2 the object over which the prediction will happend
#' @param points the data points over which F2 will be predicted. This
#' should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.F2 <- function(f2, points, spectra) {
  # get the fns for the given spectrums
  # spectra is a list of type list(list(t = 0, spectra = list(s, s, s,  list(t = 0, spectra = list(s, s, s)))
  # we need to take only the t=0 one
  cat('js = ', unlist(lapply(spectra[[1]]$spectra[[1]], function(spec) spec$js)))
  spectraForTZero <- Filter(function(s) s$t == 0, spectra)[[1]]$spectra
  fns <- getFns(f2, points, spectraForTZero)
  fns <- fns[-match(c('Q2', 'x'), names(points))]  # remove the Q2 and x columns if any
  # fit the linear model
  f2Data <- expVal(f2)
  fit <- lm(f2Data ~ ., data = fns)
  # finally return the predicted value for this linear fit for the
  # points used in the fit
  val <- predict(fit)
  if(is.null(attr(f2, 'complete')))
    val
  else
    list(val = val, coefficients = summary(fit)$coefficients)
}
#' This function returns a term in the sum for F2 with the
#' exception of the overall constant.
#' @export
#' @return A data frame where the first columns are the points passed
#' and the next ones are the values of fns for each one of the kernels
getFns.F2 <- function(f2, points, spectraForTZero) {
  newColumns <- data.frame(Q2 = points$Q2, x = points$x)
  # iterate over each kernel's spectrum
  lapply(spectraForTZero, function(s) {
    # iterate over each Reggeon for the given spectrum
    # remember, the tr1 and tr2 are not data about the reggeons
    lapply(s[names(s) == ''], function(spec) {
      fn <- apply(points, 1, function(row) {
        Q2 <- row[1]
        x  <- row[2]
        fN(Q2, x, spec$js, spec$wf)
      })
      newColumns <<- cbind(newColumns, fn)
    })
  })
  newColumns
}

fN <- Vectorize(function(Q2, x, J, wf) {
  As <- get('As', envir = ihqcdEnv)
  z  <- get('z', envir = ihqcdEnv)
  h  <- get('h', envir = ihqcdEnv)
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getMode(Q2 = Q2, h = h)$factor
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  x^(1 - J) * Q2^J * integral$value
}, c('Q2', 'x'))

#' @export
expVal.F2 <- function(f2) f2$data$F2
#' @export
expErr.F2 <- function(f2) f2$data$err
#' @export
expKinematics.F2 <- function(f2) as.data.frame(f2$data[c('Q2', 'x')])
