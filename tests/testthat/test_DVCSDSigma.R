context('DVCSDSigma')

init()

test_that('Experimental data is being loaded', {
  dvcs <- DVCSDSigma()
  expect_equal(length(expVal(dvcs)), 76)
  expect_equal(length(expErr(dvcs)), 76)
  expect_equal(names(expKinematics(dvcs)),c('Q2', 'W', 't'))
})

test_that('Computation of the differential cross-section is working.', {
  k     <- kernelUnit(UJgTest, kernelName = 'g')
  s     <- k$findKernel()
  km    <- kernelUnit(UJgTest, kernelName = 'm')
  sm    <- km$findKernel()
  gs   <- list(
            list(list(g0 = 1.2, g1 = 1), list(g0 = 0.2, g1 = 1), list(g0 = 0.3, g1 = 1)),
            list(list(g0 = 0.1, g1 = 1), list(g0 = 0.1, g1 = 1), list(g0 = 10, g1 = 1))
          )
  dvcs <- DVCSDSigma()
  # some dummy spectra
  spectra <- lapply(getNeededTVals(dvcs), function(t) list(t = t, spectra = list(s, sm)))
  p <- predict(dvcs, gs = gs, spectra = spectra)
  expect_that(p,  is_a('numeric'))
  expect_equal(length(p),  length(loadData(dvcs)$Q2))
})

test_that('Izs are being computed', {
  dvcs  <- DVCSDSigma()
  k     <- kernelUnit(UJgTest, kernelName = 'g')
  s     <- k$findKernel()
  km    <- kernelUnit(UJgTest, kernelName = 'm')
  sm    <- km$findKernel()
  # some dummy spectra
  spectra <- lapply(getNeededTVals(dvcs), function(t) list(t = t, spectra = list(s, sm)))
  Izs <- getIzs(dvcs, spectra = spectra)
})

test_that('Plot works', {
  k     <- kernelUnit(UJgTest, kernelName = 'g', numReg = 4)
  dvcsds  <- DVCSDSigma()
  points <- rbind(enlargeDataWithTs(dvcsds, ts = seq(-1, 0, len = 31)),
                  expKinematics(dvcsds))
  # spectra
  spectra <- lapply(unique(points$t),
                    function(t)
                      list(t = t,
                           spectra = list(k$findKernel(
                              .t = t,
                              invls = 5.707847,
                              a = -4.368889,
                              b = 0.614573,
                              c = 0.765144,                                                                               d = 0.000614))))
  # get the Izs for the experimental values
  expIzs <- getIzs(dvcsds, spectra = spectra)
  # find the best gs
  gs <- getBestGs(dvcsds, Izs = expIzs, numGs = 8)
  Izs <- getIzs(dvcsds, spectra = spectra, points = points)
  pred <- cbind(points,
                predicted = predict(dvcsds, Izs = Izs, points = points, gs = gs))
  plot(dvcsds, predicted = pred, numGraphs = 4)
})

