context('IHQCD')

init()

test_that('The ihqcdEnv is being properly injected by the function injectIHQCDEnv', {
  # if we just do environment(f) <- ihqcdEnv then we won't have access to possible defined
  # variables, we need to merge the environments.
  # this test function has a tt variable which is not part of the global environment,
  # so tt can't be accessed from ihqcdEnv only
  t1 <- function() {
    tt <- 10
    test <- injectIHQCDEnv(function() {
      tt * A0
    })
    test()
  }
  expect_equal(t1(), 10 * get('A0', envir = ihqcdEnv))
})
