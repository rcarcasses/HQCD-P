context('KernelUnit')

# this setup the ihqcdEnv, required for the tests
init()

test_that('kernelUnit can be computed for different values of t', {
  # create the kernel
  k <- kernelUnit(UJgTest)
  #k$setNewPotential(UJgTest, c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  kVal <- k$findKernel(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117)
  # make sure the intercepts are the right one for t=0
  expect_lt(abs(kVal[[1]]$js - 1.17), 1e-2)
  expect_lt(abs(kVal[[2]]$js - 1.08), 1e-2)
  expect_lt(abs(kVal[[3]]$js - 0.97), 1e-2)
  # make sure the intercepts are the right one for t=2.5
  kVal <- k$findKernel(.t = 2.5, invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117)
  expect_lt(abs(kVal[[1]]$js - 1.58), 1e-2)
  expect_lt(abs(kVal[[2]]$js - 1.26), 1e-2)
  expect_lt(abs(kVal[[3]]$js - 1.09), 1e-2)
})
