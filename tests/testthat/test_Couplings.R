context('Couplings')

init()


test_that('Coupling functions are being loaded and computing the expected results', {
  expect_equal(kU1(1,1), 128 * sqrt(pi)/3)     # The numerical value was computed in Wolfram Mathematica
  expect_equal(kPhi(1,1), 128 * sqrt(pi))       # The numerical value was computed in Wolfram Mathematica
})
