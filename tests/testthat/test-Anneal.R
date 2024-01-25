# anneal
test_that("anneal iterates k's", {
  data = tibble(
    x = 1:3,
    y = 1:3,
    datetime = as.Date("2017-01-01") + 0:2,
  ) %>%
    as_tsibble()

  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, new_data = data$x)
  data = data %>% mutate(fitted_obs = pred_lm)
  d = digest(data, 1, 1)

  out = anneal(
    data = data,
    fitted_obs = "fitted_obs",
    digest = d,
    resolution = 1,
    range_start = 0,
    range_end = 0,
    loess_fit = fit_lm,
    loss_fn = rmse
  )

  expected_fragments = tibble(
    idx_upsam = c(c(2,3),c(1,2,3)),
    adj_idx = c(c(3,4),c(3,4,5)),
    k = c(c(1,1),c(2,2,2)),
    .pred_obs = setNames(idx_upsam, c(c(1,2),c(1,2,3))),  # ref: https://adv-r.hadley.nz/vectors-chap.html#attr-names
    final_idx = adj_idx,
    shift = 0,
    datetime = c(
      as.Date("2017-01-03") + 0:1,
      as.Date("2017-01-03") + 0:2
    )
  )
  expected_losses = tibble(
    k = c(1,2),
    shift = c(0,0),
    loss = c(1,2)
  )

  expect_equal(out$losses, expected_losses)
  expect_equal(out$fragments, expected_fragments)
})

test_that("anneal", {
  data = tibble(
    x = 1:8,
    y = 1:8,
    datetime = as.Date("2017-01-01") + 0:7,
    # idx_data = 1:8,
    ) %>%
    as_tsibble()

  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, new_data = data$x)
  data = data %>% mutate(fitted_obs = pred_lm)
  d = digest(data, 2, 4)

  expected_fragments = tibble(
    idx_upsam = rep(3:8, 3),
    adj_idx = idx_upsam + 4,
    k = 1,
    .pred_obs = setNames(idx_upsam, rep(1:6,3)),  # ref: https://adv-r.hadley.nz/vectors-chap.html#attr-names
    final_idx = c(6:11,7:12,8:13),
    shift = c(rep(-1,6),rep(0,6),rep(1,6)),
    datetime = c(
      as.Date("2017-01-06") + 0:5,
      as.Date("2017-01-07") + 0:5,
      as.Date("2017-01-08") + 0:5
    )
  )
  expected_losses = tibble(
    k = c(1,1,1),
    shift = c(-1,0,1),
    loss = c(3,4,5)
  )

  out = anneal(
    data = data,
    fitted_obs = "fitted_obs",
    digest = d,
    resolution = 1,
    range_start = -1,
    range_end = 1,
    loess_fit = fit_lm,
    loss_fn = rmse
    )
  expect_equal(out$losses, expected_losses)
  expect_equal(out$fragments, expected_fragments)
})


# digest
test_that("digest handles tsibbles", {
  data = tibble(
      x = 101:112,
      datetime = as.Date("2017-01-01") + 0:11
    ) %>%
    as_tsibble()
  expected = tibble(
    x = c(107:112, 103:112),
    datetime = c(as.Date("2017-01-07") + 0:5, as.Date("2017-01-03") + 0:9),
    idx = c(7:12, 3:12),
    k = c(rep(1,6), rep(2,10)),
    adj_idx = c(11:16, 11:20)
  )
  out = digest(data, 2, 4)
  expect_equal(out, expected)
})

test_that("digest, no partial overlaps", {
  data = tibble(x = 101:112)
  expected = tibble(
    x = c(107:112, 103:112),
    idx = c(7:12, 3:12),
    k = c(rep(1,6), rep(2,10)),
    adj_idx = c(11:16, 11:20)
  )
  out = digest(data, 2, 4)
  expect_equal(out, expected)
})

test_that("digest, partial overlaps included", {
  data = tibble(x = 101:112)
  expected = tibble(
    x = c(104:112, 101:112),
    idx = c(4:12, 1:12),
    k = c(rep(1,9), rep(2,12)),
    adj_idx = c(8:16, 9:20)
  )
  out = digest(data, 5, 4)
  expect_equal(out, expected)
})

test_that("digest, partial overlaps excluded", {
  data = tibble(x = 101:112)
  expected = tibble(
    x = c(104:112),
    idx = c(4:12),
    k = c(rep(1,9)),
    adj_idx = c(8:16)
  )
  out = digest(data, 5, 4, include_partial_overlap = FALSE)
  expect_equal(out, expected)
})


test_that("digest, n_future_steps limits fragment size", {
  data = tibble(x = 101:112)
  expected = tibble(
    x = c(104:111, 101:107),
    idx = c(4:11, 1:7),
    k = c(rep(1,8), rep(2,7)),
    adj_idx = c(8:15, 9:15)
  )
  out = digest(data, 5, 4, n_future_steps = 3)
  expect_equal(out, expected)
})

test_that("fragment trimming", {
  data = tibble(
    x = c(NA,NA,2,2,NA,NA,1,4,NA),
    k = c(1,1,1,1,1,2,2,2,2)
  )
  expected_default = tibble(
    x = c(2,2,1,4),
    k = c(1,1,2,2)
  )
  expected_left_only = tibble(
    x = c(2,2,NA,1,4,NA),
    k = c(1,1,1,2,2,2)
  )
  expected_right_only = tibble(
    x = c(NA,NA,2,2,NA,1,4),
    k = c(1,1,1,1,2,2,2)
  )
  expected_neither = data

  expect_equal(data %>% trim_fragments_na(x), expected_default)
  expect_equal(data %>% trim_fragments_na(x, right = FALSE), expected_left_only)
  expect_equal(data %>% trim_fragments_na(x, left = FALSE), expected_right_only)
  expect_equal(data %>% trim_fragments_na(x, left = FALSE, right = FALSE), expected_neither)
})

# calc_ortho_vec
test_that("flat ts's ortho vec points up", {
  data = tibble(x = c(0,0,0,0,0))
  expected = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 4)
  expect_equal(out, expected)
})

test_that("isosceles ts's ortho vec points up", {
  data = tibble(x = c(0,5,0))
  expected = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 2)
  expect_equal(out, expected)
})

test_that("y=x ts's ortho vec points at y=(-x)", {
  data = tibble(x = c(0,1,2))
  expected = tibble(x = cos(3/4*pi), y = sin(3/4*pi))
  out = calc_ortho_vec(data$x, 2)
  expect_equal(out, expected)
})

test_that("weighting ortho vec works", {
  data = tibble(x = c(0,1,2,1))
  expected = tibble(x = 0, y = 1)
  out = calc_ortho_vec(data$x, 2, wts = c(1,1,2))
  expect_equal(out, expected)
})

test_that("y=(1/2)x w/ scaler = 2 ts's ortho vec points at y=(-x)", {
  data = tibble(x = c(0,0.5,1))
  expected = tibble(x = cos(3/4*pi), y = sin(3/4*pi))
  out = calc_ortho_vec(data$x, 2, scaler = 2)
  expect_equal(out, expected)
})


test_that("weighed rmse", {
  a = c(1,1,1,5,2)
  b = c(1,3,3,3,3)
  expect_equal(weighted_rmse(a, b), 2)
})
