test_that("rmse", {
  a = c(1,1,3,5)
  b = c(1,3,1,1)
  expect_equal(rmse(a, b), sqrt(6))
})

test_that("slope-weighed rmse", {
  a = c(1,1,1,5,2)
  b = c(1,3,3,3,3)
  expect_equal(slope_weighted_rmse(a, b), 2)
})

test_that("recency-weighed rmse", {
  a = c(1,1,1)
  b = c(1,5,3)
  # so weights should be c(0,0.5,1)
  expect_equal(recency_weighted_rmse(a, b), 2)
})
