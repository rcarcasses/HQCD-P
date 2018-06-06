DSigma <- function(procName, postName = '') {
  # add the generic DSigma class
  obs <- ProcessObservable(paste0(procName, 'DSigma', postName))
  class(obs) <- append(class(obs), 'DSigma', after = length(class(obs)) - 1)
  obs
}

# Conversion factor between GeV^2 and mili/nano barn units
GEVMinus2ToNB <- 3.894*10^5
GEVMinus2ToMB <- 0.3894

#' @export
getNeededTVals.DSigma <- function(x) unique(expKinematics(x)$t)

getAmplitude <- function(x, ...) UseMethod('getAmplitude')
getAmplitude.DSigma <- function(dsigma, Izs, IzsBar, points, ...) {
  # if it is a composed object then just call getAmplitude on each children
  if(!is.null(dsigma$tt) && !is.null(dsigma$ll))
    return(list(tt = getAmplitude(dsigma$tt, Izs, IzsBar, points = points, ...),
                ll = getAmplitude(dsigma$ll, Izs, IzsBar, points = points, ...)))

  rowSums(Izs * IzsBar, na.rm = TRUE)
}

getAmplitudeNMC1 <- function(x, ...) UseMethod('getAmplitudeNMC1')
getAmplitudeNMC1.DSigma <- function(dsigma, Izs, IzsBar, points, ...){if(!is.null(dsigma$tt) && !is.null(dsigma$ll))
  return(list(tt = getAmplitudeNMC1(dsigma$tt, Izs, IzsBar, points = points, ...),
              ll = getAmplitudeNMC1(dsigma$ll, Izs, IzsBar, points = points, ...)))
  rowSums(Izs * IzsBar, na.rm = TRUE)
  }

getAmplitudeNMC2 <- function(x, ...) UseMethod('getAmplitudeNMC2')
getAmplitudeNMC2.DSigma <- function(dsigma, Izs, IzsBar, points, ...){
  if(!is.null(dsigma$tt) && !is.null(dsigma$ll))
    return(list(tt = getAmplitudeNMC2(dsigma$tt, Izs, IzsBar, points = points, ...),
                ll = getAmplitudeNMC2(dsigma$ll, Izs, IzsBar, points = points, ...)))
  rowSums(Izs * IzsBar, na.rm = TRUE)
}

#' Get the integrals on z
#' @export
getIzs.DSigma <- function(dsigma, spectra, points) {
  # if it is a composed object then just call getIzs on each children
  if(!is.null(dsigma$tt) && !is.null(dsigma$ll))
    return(list(tt = getIzs(dsigma$tt, spectra, points = points),
                ll = getIzs(dsigma$ll, spectra, points = points)))
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s, function(spec) paste0('fn.', spec$name)))))
  pb <- progress_bar$new(format = " getIzs.DSigma [:bar] :percent eta: :eta",
                          total = length(points[[1]]), clear = FALSE, width= 60)
  df <- as.data.frame(Reduce(rbind, mclapply(apply(points, 1, as.list), function(row) {
    pb$tick()
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
        # use the definition passed with the object to compute this integral
        # if any, the attribute has to be a function with the right arity
        if(!is.null(attr(dsigma, 'IzN')))
          attr(dsigma, 'IzN')(dsigma, W, Q2, spec)
        else
          IzN(dsigma, W, Q2, spec$js, spec$wf)
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  }, mc.cores = cores)))
  attr(df, 'row.names') <- 1:length(df[[1]])
  df
}

#' Get the integrals on z bar
#' @export
getIzsBar.DSigma <- function(dsigma, spectra, points, zstar, hpars) {
  # if it is a composed object then just call getIzs on each children
  if(!is.null(dsigma$tt) && !is.null(dsigma$ll))
    return(list(tt = getIzs(dsigma$tt, spectra, points = points),
                ll = getIzs(dsigma$ll, spectra, points = points)))
  fnNames <- unlist(lapply(spectra[[1]]$spectra,
                           function(s)
                             unlist(lapply(s, function(spec) paste0('fn.', spec$name)))))
  # TODO improve this, right now I'm just caching the computation
  # maybe we should rewrite all the code to make it more clear and efficient
  pb <- progress_bar$new(format = " getIzsBar.DSigma [:bar] :percent eta: :eta",
                          total = length(points[[1]]), clear = FALSE, width= 60)
  IzBarCache <- list()
  df <- as.data.frame(Reduce(rbind, lapply(apply(points, 1, as.list), function(row) {
    pb$tick()
    row <- as.list(row)
    t  <- row$t
    # get the spectra of all kernels for a given value of t
    spectraForT <- Filter(function(s) s$t == t, spectra)[[1]]$spectra
    # iterate over each kernel's spectrum
    r <- unlist(lapply(spectraForT, function(s) {
      # s: spectrum of a single kernel, have many reggeons
      # iterate over each Reggeon for the given spectrum
      lapply(s, function(spec) {
        # here we use the fact that the integral value depends only on
        # t and the reggeon in case therefore we cache the result
        id <- paste0(t, spec$index)
        if(is.null(IzBarCache[[id]]))
          IzBarCache[[id]] <- if(!is.null(attr(dsigma, 'IzNBar')))
            # use the definition passed with the object to compute this integral
            # if any, the attribute has to be a function with the right arity
            attr(dsigma, 'IzNBar')(dsigma, W, Q2, spec, zstar, hpars)
          else
            IzNBar(dsigma, spec$js, spec$wf, spec$dJdt, zstar, hpars)

        IzBarCache[[id]]
      })
    }), recursive = TRUE)
    names(r) <- fnNames
    r
  })))
  attr(df, 'row.names') <- 1:length(df[[1]])
  df
}

#' @export
IzN.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalStateFactor(dsigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  W^(2*J) *integral$value
}

#' @export
IzNNMC1.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J + 1.5) * As))
  t2fun <- getExternalStateFactor(dsigma, Q2 = Q2)
  t3fun <- splinefun(z, DperpPsi(wf))
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  W^(2*J) *integral$value
}

#' @export
IzNNMC2.DSigma <- function(dsigma, W, Q2, J, wf) {
  t1fun <- splinefun(z, exp((-J - 0.5) * As))
  t2fun <- getExternalStateFactor(dsigma, Q2 = Q2)
  t3fun <- splinefun(wf$x, wf$y)
  integral <- integrate(function(x)  t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
  # return the full thing needed for the amplitude
  -W^(2*J) *integral$value
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
enlargeData.DSigma <- function(dsigma, ts = seq(-1, 0, 0.05)) {
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
