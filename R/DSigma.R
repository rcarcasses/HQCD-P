DSigma <- function(procName) {
  # add the generic DSigma class
  obs <- ProcessObservable(paste0(procName, 'DSigma'))
  class(obs) <- append(class(obs), 'DSigma', after = length(class(obs)) - 1)
  obs
}

# Conversion factor between GeV^
GEVMinus2ToNB <- 3.894*10^5

#' @export
getNeededTVals.DSigma <- function(x) unique(expKinematics(x)$t)

getAmplitude <- function(x, ...) UseMethod('getAmplitude')
getAmplitude.DSigma <- function(dsigma, fns, gs, points, ...) {
  # compute g(t) for each corresponding fn, this return a dataframe
  # were each column is the value (a vector) of g(t) for the given values
  # of t and the n column is the g(t) of the n fn in the fns dataframe
  # see the test file for some explanation 'tests/testthat/test_DSigma.R'
  gts <- apply(gs, 1, function(row) {
    rowSums(t(row * t(outer(points$t, 0:(length(gs) - 1), `^`))))
  })
  rowSums(fns * gts, na.rm = TRUE)
}

#' Get fns times dJdt
#' @export
getFns.DSigma <- function(dsigma, spectra, points) {
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s, function(spec) paste0('fn.', spec$name)))))
  df <- data.frame(row.names = fnNames)
  as.data.frame(t(apply(points, 1, function(row) {
    row <- as.list(row)
    Q2 <- row$Q2
    W  <- row$W
    t  <- row$t
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    r <- unlist(lapply(spectraForT, function(s) {
      # s: spectrum of a single kernel, have many reggeons
      # iterate over each Reggeon for the given spectrum
      lapply(s, function(spec) {
        fN(dsigma, W, Q2, spec$js, spec$wf)
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  })))
}

fN.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalStateFactor(dsigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  (1 - 1i/ tan(pi * J / 2)) * W^(2*J) *integral$value
}

#' @export
getBestGs.DSigma <- function(dsigma, fns, numGs, startGs = NULL) {
  # first we need to define an function depending only of the gs
  # to be optimized
  fn <- function(gs) rss(dsigma, fns = fns, gs = as.data.frame(matrix(gs, ncol = numGs / length(fns))))
  # if is the fist time just put something there
  if(is.null(startGs))
    startGs <- rep(1, len = numGs)

  op <- optim(startGs,
              fn = fn, method = 'BFGS', hessian = FALSE, control = list(maxit = 1000))
  # store the best gs found so they can be used as a starting point of the next call
  flog.debug(paste('DSigma bestGs  =', do.call(paste, as.list(format(op$par, digits = 4))), ' in', op$counts[1], ' steps'))
  as.data.frame(matrix(op$par, ncol = numGs / length(fns)))
}

getExternalStateFactor <- function(x, ...) UseMethod('getExternalStateFactor')

getExternalStateFactor.default <- function(x, ...) 'getExternalStateFactor have to be implemented for this process'

#' @export
expVal.DSigma <- function(dsigma) dsigma$data$dsigma
#' @export
expErr.DSigma <- function(dsigma) dsigma$data$deltaDSigma
#' @export
expKinematics.DSigma <- function(dsigma) dsigma$data[c('Q2', 'W', 't')]

#' For these processes the enlarge data introduces many values of t for
#' each one of the Q2 and W combinations available
#' @export
enlargeData.DSigma <- function(dsigma, ts = seq(-1, 0, len = 11)) {
  # get all the different combinations of Q2s and Ws
  QsAndWs <- unique(dsigma$data[,c('Q2', 'W')])
  as.data.frame(lapply(as.data.frame(Reduce(function(acc, r) {
    val <- rbind(acc, r)
    rownames(val) <- NULL
    val
  },unlist(
    apply(QsAndWs, 1,
          function(row)
            lapply(ts, function(t) c(list(t = t), row))),
    recursive = FALSE)
  )), `mode<-`, 'numeric'))
}
