context('VMPDSigma')

init()

test_that('Experimental data is being loaded', {
  vmp <- VMPDSigma()
  expect_equal(length(expVal(vmp)), 120)
  expect_equal(length(expErr(vmp)), 120)
  expect_equal(names(expKinematics(vmp)),c('m', 'Q2', 'W', 't'))
})
