context('ppDSigma')

init()

test_that('Experimental data is being loaded', {
  pp <- ppDSigma()
  expect_equal(length(expVal(pp)), 650)
  expect_equal(length(expErr(pp)), 650)
  expect_equal(names(expKinematics(pp)),c('Q2', 'W', 't'))
})
