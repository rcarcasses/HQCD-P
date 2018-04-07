context('VMPDSigma')

init()

test_that('Experimental data is being loaded', {
  vmp <- VMPDSigma('Rho0')
  expect_equal(length(expVal(vmp)), 35)
  expect_equal(length(expErr(vmp)), 35)
  expect_equal(names(expKinematics(vmp)),c('Q2', 'W', 't'))
})
