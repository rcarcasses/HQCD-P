context('HQCDP')

init()
clearCache()

test_that('HQCDP addKernel and addProcessObservable are functional', {
  # test the addKernel
  p <- HQCDP()
  p <- addKernel(p, potential = UJgTest,
                  numReg = 2,
                  comment = 'Leading twist gluon sector',
                  kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                  optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  p <- addKernel(p, potential = UJgTest,
                    numReg = 2,
                    comment = 'Leading twist gluon sector',
                    kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                    optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  expect_true(length(p$kernels) == 2)
})

test_that('HQCDP fit works', {
  # test the addKernel
  p <- HQCDP()
  p <- addKernel(p, potential = UJgTest,
                 numReg = 2,
                 comment = 'Leading twist gluon sector',
                 kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                 optimPars = c(invls = 1/0.153, a = -4.35))
  p <- addKernel(p, potential = UJgTest,
                 numReg = 2,
                 comment = 'Leading twist gluon sector',
                 kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                 optimPars = c(invls = 1/0.153, c = 0.626, d = -0.117))
  expect_true(length(p$kernels) == 2)
})

test_that('HQCDP addKernel is functional', {
  # test the addKernel
  # create a mockup f2 object
  f2   <- F2()
  dvcs <- DVCSDSigma()
  p <- HQCDP()
  p <- addProcessObservable(p, f2)
  p <- addProcessObservable(p, dvcs)
  expect_true(length(p$processes) == 2)
})

#test_that('The same result as in https://arxiv.org/abs/1704.08280 is gotten', {
#  # compute this with 400 cheb points
#  p <- HQCDP()
#  p <- addKernel(p)
#  p <- addProcessObservable(p, F2())
#  expect_equal(rss(p, pars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117)), - 405.7502, tolerance = 1e-3)
#})

test_that('The gs are converted from a vector to a data.frame back and forth properly', {
  # create to sample vectors of the same size
  g0 <- 1:6
  g1 <- 7:12
  # build a dataframe with them
  gs <- data.frame(g0, g1)
  # get a vector from this
  allGs <- c(g0, g1)
  # now convert it back to a dataframe
  df <- as.data.frame(matrix(allGs, ncol = 2))
  # check that df and gs are the same
  expect_true(Reduce(and, df == gs))
})

test_that('The computations of each one of the Reggeon\'s data is assembled properly', {
  # create some dummy list representing a single kernel spectrum
  r1 <- list(js = 1, wf = 3)
  r2 <- list(js = 2, wf = 6)
  r3 <- list(js = 3, wf = 2)
  r4 <- list(js = 4, wf = 4)
  r5 <- list(js = 5, wf = 5)
  r6 <- list(js = 6, wf = 7)
  s1 <- list(r1, r2)
  s2 <- list(r3, r4, r5)
  s3 <- list(r6, r1)
  s4 <- list(r2, r3, r4)
  # the spectra have to be a list of results returned by k$findKernel()
  # with as many entries as kernels are considered, here we mimic 2
  spectra <- list(list(t = 0, spectra = list(s1, s2)), list(t = 1, spectra = list(s3, s4)))
  # now suppose that we compute this using the findReggeonData function, desirable
  # for parallel computing, we need to assemble the results with the same structure
  rawSpectra  <- list(r1, r2, r3, r4, r5,
                      r6, r1, r2, r3, r4)
  # we need to pack the rawSpectra into pieces of length numReg
  # which are the amount of reggeons of each kernel
  # notice that in this example the first kernel have two reggeon, the second 3
  # and the last 1
  ts <- c(0, 1)
  numRegs <- c(2, 3)
  builtSpectra <- convertRawSpectra(rawSpectra, numRegs, ts)
  expect_true(identical(builtSpectra, spectra))
})

test_that('completeWithFixedVal works properly', {
  vals <- c(invls = 2, a = 1, b = NA, c = NA, d = 1)
  HQCDP(fixed = list(fixedPars = c(NA, NA, 1.1, 2.2, NA), zstar = 2.3)) -> p
  valsCompleted <- completeWithFixedVal(p, vals, 'fixedPars')
  expect_equal(valsCompleted,  c(invls = 2, a = 1, b = 1.1, c = 2.2, d = 1))
  # notice that we need to pass NA while replacing a number
  expect_equal(attr(p, 'fixed')$zstar,  completeWithFixedVal(p, NA, 'zstar'))
})

test_that('spectra interpolation works', {
  flog.threshold(DEBUG)
  ts <- seq(-5, 0, len = 40)
  fp <- structure(list(), class = c('FakeProcess', 'ProcessObservable'))
  # let's first compute some spectra exactly
  getNeededTVals.FakeProcess <- function(...) ts
  HQCDP() %>% addKernel %>% addProcessObservable(fp)-> p
  spectra <- getSpectra(p)
  # now let's compute the spectra using some interpolation instead
  attr(p, 'useTVals') <- seq(-5, 0, len = 10)
  interpolatedSpectra <- getSpectra(p)
  # now let's compare the results obtained by interpolation against the true results
  flattenSpectra <- unlist(lapply(spectra, `[[`, 'spectra'), recursive = FALSE)
  flattenSpectraInter <- unlist(lapply(interpolatedSpectra, `[[`, 'spectra'), recursive = FALSE)
  # compare the two 4th trajectories
  index <- 4
  Jt <- unlist(lapply(lapply(flattenSpectra, `[[`, index), `[[`, 'js'))
  JtInter <- unlist(lapply(lapply(flattenSpectraInter, `[[`, index), `[[`, 'js'))
  dJdt <- unlist(lapply(lapply(flattenSpectra, `[[`, index), `[[`, 'dJdt'))
  dJdtInter <- unlist(lapply(lapply(flattenSpectraInter, `[[`, index), `[[`, 'dJdt'))
  expect_equal(Jt, JtInter, tolerance = 0.001)
  expect_equal(dJdt, dJdtInter, tolerance = 0.001)
  # compare the 27 3th wavefunctions
  wf1 <- spectra[[27]]$spectra[[1]][[3]]$wf
  wf2 <- interpolatedSpectra[[27]]$spectra[[1]][[3]]$wf
  expect_equal(wf1$y, wf2$y, tolerance = 0.001)
  # clean the mess
  rm(getNeededTVals.FakeProcess)
})
