context('Loader')
# create some dummy vector
Avec <- (1:100)^2

test_that('Entities can be loaded from RDS files', {
  # save it as rds
  saveRDS(Avec, 'Avec.rds')
  Avec2 <- loadExt('Avec')
  expect_equal(Avec, Avec2)
  # clean up
  if (file.exists('Avec.rds'))
    file.remove('Avec.rds')
})

test_that('Entities can be loaded from CSV files', {
  # save it as csv
  write.csv(Avec, 'Avec.csv', row.names = FALSE)
  Avec3 <- loadExt('Avec', afterLoad = function(x) x[[1]])
  # notice that while reading from a CSV we get a dataframe
  expect_equal(Avec, Avec3)
  # clean up
  if (file.exists('Avec.csv'))
    file.remove('Avec.csv')
})
