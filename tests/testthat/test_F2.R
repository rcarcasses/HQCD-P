context('F2 process')

flogLevel <- flog.threshold(INFO)
init()

test_that('Izs are computed', {
  k   <- kernelUnit(UJgTest, kernelName = 'g')
  s   <- k$findKernel()
  f2  <- F2()
  # pay attention to the structure of the spectra object
  spectra <- list(list(t = 1, spectra = list(s)), list(t = 0, spectra = list(s, s, s)))
  Izs <- getIzs(f2, spectra = spectra)
  expect_equal(length(Izs), 9)
})

test_that('predict and rss generic functions are called properly', {
  # recall that predict is used by rssComplete, if this work then the former as well
  k   <- kernelUnit(UJgTest)
  s   <- k$findKernel()
  f2  <- F2()
  # pay attention to the structure of the spectra object
  spectra <- list(list(t = 1, spectra = list(s)), list(t = 0, spectra = list(s, s, s)))
  val <- rss(f2, spectra = spectra)
  # first the simple answer
  #expect_lt(abs(val - 990.0943), 1e-3)
  expect_equal(val, 1214.598, tolerance = 1e-3)
  # now let's check passing the 'complete' attribute
  attr(f2, 'complete') <- TRUE
  val <- rss(f2, spectra = spectra)
  # first the simple answer
  expect_equal(val$val, 1214.598, tolerance = 1e-3)
})

test_that('getNeededTVals works properly', {
  f2  <- F2()
  expect_true(getNeededTVals(f2) == 0)
})
