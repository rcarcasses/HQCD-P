context('DVCSSigma')

init()

test_that('Experimental data is being loaded', {
  dvcs <- DVCSSigma()
  expect_equal(length(expVal(dvcs)), 58)
  expect_equal(length(expErr(dvcs)), 58)
  expect_equal(names(expKinematics(dvcs)),c('Q2', 'W'))
  expect_equal(attr(dvcs, 'tmin'), -1.0)
  expect_equal(attr(dvcs, 'tmax'), 0.0)
})

test_that('Total cross-section can be computed', {
  HQCDP() %>%
  addKernel(potential = UJgTest,
              numReg = 2,
              comment = 'Leading twist gluon sector',
              kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
              optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117)
            ) %>%
  addProcessObservable(DVCSSigma()) -> p
  rss(p)
})
