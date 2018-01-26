context('VMPSigma')

init()

test_that('Experimental data is being loaded', {
  vmp <- VMPSigma()
  expect_equal(length(expVal(vmp)), 113)
  expect_equal(length(expErr(vmp)), 113)
  expect_equal(names(expKinematics(vmp)),c('m', 'Q2', 'W'))
})
