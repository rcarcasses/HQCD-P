context('HQCDP')

test_that('HQCDP addKernel and addProcessObservable are functional', {
  # test the addKernel
  # create a mockup f2 object
  f2 <- list(b = 10, c = 1)
  class(f2) <- append(class(f2), 'kernelUnit')
  hqcd <- HQCDP()
  hqcd <- addKernel(hqcd, f2)
  hqcd <- addKernel(hqcd, f2)
  hqcd <- addKernel(hqcd, f2)
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
