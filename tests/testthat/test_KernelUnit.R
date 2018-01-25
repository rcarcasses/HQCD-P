context('Fitter')

# Some test function (actually the one being used in the NMC paper)
UJg <- function(J, invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117) {
  # the deformation of the potential of the graviton
  u2 + (J - 2) * ((2 * e2As * invls^2) * (1 + d / l1_2) + e2A * (J + 2)
                  + (a * aF + b * bF + c * cF))
}

test_that('kernelUnit can be computed for different values of t', {
  # create the kernel
  k <- kernelUnit()
  k$setNewPotential(UJg, c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  kVal <- k$findKernel()
  # make sure the intercepts are the right one for t=0
  expect_less_than(abs(kVal[[1]]$js - 1.17), 1e-2)
  expect_less_than(abs(kVal[[2]]$js - 1.08), 1e-2)
  expect_less_than(abs(kVal[[3]]$js - 0.97), 1e-2)
  # make sure the intercepts are the right one for t=2.5
  kVal <- k$findKernel(atT = 2.5)
  expect_less_than(abs(kVal[[1]]$js - 1.58), 1e-2)
  expect_less_than(abs(kVal[[2]]$js - 1.26), 1e-2)
  expect_less_than(abs(kVal[[3]]$js - 1.09), 1e-2)
})
