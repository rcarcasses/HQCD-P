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

test_that('The same result as in https://arxiv.org/abs/1704.08280 is gotten', {
  # compute this with 400 cheb points
  p <- HQCDP()
  p <- addKernel(p, potential = UJgTest,
                 numReg = 4,
                 comment = 'Leading twist gluon sector',
                 kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                 optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  p <- addProcessObservable(p, F2())
  expect_equal(rss(p) - 405.7502, tolerance = 1e-3)
})

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
  s3 <- list(r6)
  # the spectra have to be a list of results returned by k$findKernel()
  # with as many entries as kernels are considered, here we mimic 2
  spectra <- list(s1, s2, s3)
  # now suppose that we compute this using the findReggeonData function, desirable
  # for parallel computing, we need to assemble the results with the same structure
  rawSpectra  <- list(r1, r2, r3, r4, r5, r6)
  # we need to pack the rawSpectra into pieces of length numReg
  # which are the amount of reggeons of each kernel
  # notice that in this example the first kernel have two reggeon, the second 3
  # and the last 1
  numRegs <- c(2, 3, 1)
  builtSpectra <- convertRawSpectra(rawSpectra, numRegs)
  expect_true(identical(builtSpectra, spectra))
})
