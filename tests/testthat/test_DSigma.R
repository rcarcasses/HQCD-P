context('DSigma')

test_that('gs data frame is combined with values of t properly', {
  # let's create some dummy gs dataframe
  # in this case represents the first 3 coefficients
  # so g(t) = g0 + g1 * t + g2 * t^2
  gs <- data.frame(g0 = 1:4, g1 = 5:8, g2 = 10:13)
  # create some dummy values of t
  t <- 1:10
  gts <- apply(gs, 1, function(row) {
    rowSums(t(row * t(outer(t, 0:(length(gs) - 1), `^`))))
  })

  expect_equal(sum(rowSums(gts)), 19240)
  # the outer product is very handy, here is an example of output
  # > outer(t, 0:(length(gs) - 1), `^`)
  # [,1] [,2] [,3]
  # [1,]    1    1    1
  # [2,]    1    2    4
  # [3,]    1    3    9
  # [4,]    1    4   16
  # [5,]    1    5   25
  # ...
  # the previous step give us a matrix where the columns are t^0, t^1 and t^2
  # now we would like to multiply each one of these columns by the respective
  # coefficient g, we can do it the following way
  # > t(row * t(outer(t, 0:(length(gs) - 1), `^`)))
  # [,1] [,2] [,3]
  # [1,]    1    2    3
  # [2,]    1    4   12
  # [3,]    1    6   27
  # [4,]    1    8   48
  # [5,]    1   10   75
  # ...
  # here row is 1:3, then we need to rowSum all these to get the fn for a given row (g0, g1, g2)
  # > rowSums(t(row * t(outer(t, 0:(length(gs) - 1), `^`))))
  # [1]   6  17  34  57  86 121 162 209 262 321
})
