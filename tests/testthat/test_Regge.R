context('Regge')

init()
r <- regge()

test_that('Trajectories intercepts are good', {
  r$setNewPotential(UJgTest, c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
  # get the trajectories
  tr <- r$get()
  # check some values
  expect_lt(abs(tr$j[100] - 6.5), 1e-3)
  expect_lt(abs(sum(tr$t[[100]] - c(19.25675, 34.78273, 48.36459, 61.30374))), 1e-3)
})
