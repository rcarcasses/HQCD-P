context('HQCDP')

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
