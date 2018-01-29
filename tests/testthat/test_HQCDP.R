context('HQCDP')

test_that('HQCDP addKernel and addProcessObservable are functional', {
  # test the addKernel
  hqcd <- HQCDP()
  hqcd <- addKernel(hqcd, potential = UJgTest,
                  numReg = 2,
                  comment = 'Leading twist gluon sector',
                  prefix = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                  optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  hqcd <- addKernel(hqcd, potential = UJgTest,
                    numReg = 2,
                    comment = 'Leading twist gluon sector',
                    kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
                    optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  expect_true(Reduce(and, unlist(lapply(hqcd$kernels, function(p) p$b)) == c(10, 10, 10)))
})

test_that('HQCDP addKernel is functional', {
  # test the addKernel
  # create a mockup f2 object
  f2 <- list(b = 10, c = 1)
  class(f2) <- append(class(f2), 'ProcessObservable')
  hqcd <- HQCDP()
  hqcd <- addProcessObservable(hqcd, f2)
  hqcd <- addProcessObservable(hqcd, f2)
  hqcd <- addProcessObservable(hqcd, f2)
  expect_true(Reduce(and, unlist(lapply(hqcd$processes, function(p) p$b)) == c(10, 10, 10)))
})
