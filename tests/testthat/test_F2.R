context('F2 process')

init()

test_that('Fns are computed', {
  k   <- kernelUnit()
  s   <- k$findKernel()
  fns <- getFns(f2, data.frame(Q2 = c(1.2), x = c(1e-4)), list(s))
  fns <- fns[c(-1,-2)]
  expect_equal(length(fns), 3)
})

test_that('predict and rssComplete generic functions are called properly', {
  # recall that predict is used by rssComplete, if this work then the former as well
  k   <- kernelUnit()
  s   <- k$findKernel()
  rss <- rssComplete(f2, spectra = list(s))
  expect_lt(abs(rss$rss - 607.2166), 1e-3)
})
