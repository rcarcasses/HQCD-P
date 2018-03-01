context('F2 process')

flogLevel <- flog.threshold(INFO)
init()

test_that('Fns are computed', {
  k   <- kernelUnit(UJgTest, kernelName = 'g')
  s   <- k$findKernel()
  f2  <- F2()
  # pay attention to the structure of the spectra object
  spectra <- list(list(t = 1, spectra = list(s)), list(t = 0, spectra = list(s, s, s)))
  fns <- getFns(f2, spectra = spectra)
  expect_equal(length(fns), 9)
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

test_that('getBestGs is working properly', {
  k <- kernelUnit(potential = UJgTest,
                  numReg = 4,
                  comment = 'Leading twist gluon sector',
                  kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                  optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117)
  )
  s <- k$findKernel()
  spectra <- list(list(t = 1, spectra = list(s)), list(t = 0, spectra = list(s)))
  f2 <- F2()
  fns <- getFns(f2, spectra = spectra)
  # find rss only specifying the spectra
  rss1 <- rss(f2, spectra)
  # find rss now with fns and gs
  rss2 <- rss(f2, fns = fns, gs = c(0.1787147, 0.1041425, 0.1853301, -1.530063))
  # now find rss using the HQCDP getBestGs
  p <- HQCDP()
  p <- addKernel(p, potential = UJgTest,
                 numReg = 4,
                 comment = 'Leading twist gluon sector',
                 kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                 optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  p <- addProcessObservable(p, f2)
  rss3 <- rss(p)
  expect_equal(rss1 * rss2 / rss3^2, 1)
})

test_that('Results of the NMC paper are obtained', {
  k <- kernelUnit(potential = UJgTest,
          numReg = 4,
          comment = 'Leading twist gluon sector, NMC case',
          kernelName = 'gluon',
          optimPars = c(invls = 5.91662420, a = -4.50176652, b = -0.63141639, c = 0.82567207, d = -0.01441428))
  s <- k$findKernel()
  spectra <- list(list(t = 1, spectra = list(s)), list(t = 0, spectra = list(s)))
  f2 <- F2()
  # add the non minimal coupling alpha value
  attr(f2, 'alpha') <- -0.0672
  fns <- getFns(f2, spectra = spectra)
  rss1 <- rss(f2, fns = fns)
  chi2 <- rss1 / (length(expVal(f2)) - 8 - 5)

})

test_that('getNeededTVals works properly', {
  f2  <- F2()
  expect_true(getNeededTVals(f2) == 0)
})
