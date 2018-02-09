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

test_that('Parallel computations can be done (Incomplete)', {
  t <- function(J, n) {
    args <- as.list(c(J = J))
    u2j  <- do.call(UJgTest, args)
    s <- computeSpectrum(z, u2j, n)
    list(t = s$energies[[n]], wf = s$wf[[n]], u2j = u2j, energies = s$energies)
  }

  .t <- 0
  kernelName <- 'test'
  tvec <- function(J, n) sapply(J, function(j) t(j, n)$t)
  Js <- seq(0.2, 2, len = 20)
  getIntercept <- function(n) {
    tspline <- function(j) tvec(j, n) - .t
    roots <- uniroot.all(tspline, c(0.2, 2.2), n = 6)
    if(length(roots) > 0)
      js <- max(roots)
    else
      js <- 0

    tresult <- t(js, n)
    list(js = js, wf = tresult$wf, u2j = tresult$u2j, name = paste0(kernelName, '.', n))
  }

  intercepts <- mclapply(1:3, getIntercept)
  expect_equal(unlist(lapply(intercepts, `[[`, 'js')), c(1.16911807113316, 1.08376877103049, 0.965954406166295))
})
