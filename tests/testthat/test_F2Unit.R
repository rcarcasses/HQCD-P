context('F2Unit')

# this setup the ihqcdEnv, required for the tests
init()

test_that('Generic functions are being invoked in F2Unit objects', {
  # create a default f2unit object
  f2unit <- F2Unit()
  # and invoke the generic function getAllFns and make sure its result is not zero
  expect_equal(length(getAllFns(f2unit)), 6)
})
