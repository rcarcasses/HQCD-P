context('DVCSDSigma')

init()

test_that('Experimental data is being loaded', {
  dvcs <- DVCSDSigma()
  expect_equal(length(expVal(dvcs)), 76)
  expect_equal(length(expErr(dvcs)), 76)
  expect_equal(names(expKinematics(dvcs)),c('Q2', 'W', 't'))
})

test_that('Computation of the differential cross-section is working.',{
  k    <- kernelUnit(UJgTest)
  s    <- k$findKernel()
  gs   <- list(
            list(list(g0 = 1, g1 = 1), list(g0 = 1, g1 = 1), list(g0 = 1, g1 = 1)),
            list(list(g0 = 1, g1 = 1), list(g0 = 1, g1 = 1), list(g0 = 1, g1 = 1))
          )
  dvcs <- DVCSDSigma()
  # some dummy spectra
  spectra <- lapply(getNeededTVals(dvcs), function(t) list(t = t, spectra = list(s, s)))
  p <- predict(dvcs, gs = gs, spectra = spectra)
  expect_that(p,  is_a('numeric'))
  expect_equal(length(p),  length(loadData(dvcs)$Q2))
})
