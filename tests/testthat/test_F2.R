context('F2 process')

init()

test_that('Fns are computed', {
  k   <- kernelUnit()
  s   <- k$findKernel()
  fns <- getFns(f2, data.frame(Q2 = c(1.2), x = c(1e-4)), list(s))
  fns <- fns[c(-1,-2)]
  expect_equal(length(fns), 3)
})
