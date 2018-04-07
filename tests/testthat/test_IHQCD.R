context('IHQCD')

init()

test_that('The IHQCD results have been globally injected', {
  expect_equal(A0, 5)
})
