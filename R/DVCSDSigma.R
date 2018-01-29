#' @export
DVCSDSigma <- function() ProcessObservable('DVCSDSigma')

#' @export
getNeededTVals.DVCSDSigma <- function(x) unique(expKinematics(x)$t)

#' Predicts the values of F2 for the points passed
#' @param dvcsDSigma the object over which the prediction will happend
#' @param points the data points over which dsigma/dt for dvcs will be predicted.
#' This should be a data frame with the same structure as the one returned
#' by expKinematics. Is important to keep the order of the columns.
#' @param spectra a collection of spectrum of different kernels which can have different amount of Reggeons, etc.
#' @export
predict.DVCSDSigma <- function(dvcsDSigma, points, spectra) {
  # get the fns for the given spectrums
  fns <- getFns(dvcsDSigma, points, spectra)
  fns <- fns[-match(c('Q2', 'W', 't'), names(points))]  # remove the Q2 and x columns if any
  # fit the linear model
  dvcsDSigmaData <- expVal(dvcsDSigma)
  fit <- lm(dvcsDSigmaData ~ ., data = fns)
  # finally return the predicted value for this linear fit for the
  # points used in the fit
  predict(fit)
}

fN <- function(Q2, W, t, J, wf) {
  As <- get('As', envir = ihqcdEnv)
  z  <- get('z', envir = ihqcdEnv)
  h  <- get('h', envir = ihqcdEnv)
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t1barfun <- splinefun(z,exp((-J + 3.5)*As))
  t2fun <- getMode(Q2 = Q2, h = h)$factor
  t22fun <- getMode(Q2 = 0, h = h)$factor
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x) t1fun(x) * t1barfun(x) * t2fun(x) * t22fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  result <- - (kU1(J)*kPhi(J)/2^J) * W^(2 * J) * (pi/2) * integral$value
}

#' @export
expVal.DVCSDSigma <- function(dvcs) dvcs$data$dsigma
#' @export
expErr.DVCSDSigma <- function(dvcs) dvcs$data$deltaDSigma
#' @export
expKinematics.DVCSDSigma <- function(dvcs) dvcs$data[c('Q2', 'W', 't')]
