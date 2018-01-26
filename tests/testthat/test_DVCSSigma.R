context('DVCSSigmaUnit')

test_that('Experimental data is being loaded', {
  dvcs <- DVCSSigmaUnit()
  expect_equal(length(expVal(dvcs)), 58)
  expect_equal(length(expErr(dvcs)), 58)
  expect_equal(names(expKinematics(dvcs)),c('Q2', 'W'))
})
