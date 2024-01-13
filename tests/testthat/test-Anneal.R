test_that("flat ts's ortho vec points up", {
  data = tibble(x = c(0,0,0,0,0))
  expc = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 4)
  expect_equal(out, expc)
})

test_that("isosceles ts's ortho vec points up", {
  data = tibble(x = c(0,5,0))
  expc = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 2)
  expect_equal(out, expc)
})

test_that("y=x ts's ortho vec points at y=(-x)", {
  data = tibble(x = c(0,1,2))
  expc = tibble(x = cos(3/4*pi), y = sin(3/4*pi))
  out = calc_ortho_vec(data$x, 2)
  expect_equal(out, expc)
})

test_that("weighting ortho vec works", {
  data = tibble(x = c(0,1,2,1))
  expc = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 2, wts = c(1,1,2))
  expect_equal(out, expc)
})
