context('DVCSDSigmaUnit')

test_that('Experimental data is being loaded', {
  dvcs <- DVCSDSigmaUnit()
  expect_equal(length(expVal(dvcs)), 76)
  expect_equal(length(expErr(dvcs)), 76)
  expect_equal(names(expKinematics(dvcs)),c('Q2', 'W', 't'))
})
