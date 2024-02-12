# anneal
test_that("anneal iterates k's", {
  data = tibble(
    x = 1:3,
    y = 1:3,
    datetime = as.Date("2017-01-01") + 0:2,
  ) %>%
    as_tsibble()

  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, newdata = data %>% select(x))
  data = data %>% mutate(fitted_obs = pred_lm)
  d = data %>% digest(
    datetime, y, 1, 1,
    min_fragment_len = 1)  # overruling the default min for this test.

  out = data %>% anneal(
    fitted_obs = fitted_obs,
    digest = d,
    range_start = 0,
    range_end = 0,
    loss_fn = rmse
  )

  expected_fragments = d %>%
    mutate(
      final_idx = adj_idx,
      final_datetime = adj_datetime,
      shift = 0
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
    ) %>%
    as_tsibble()

  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, newdata = data %>% select(x))
  data = data %>% mutate(fitted_obs = pred_lm)
  d = data %>% digest(datetime, y, 2, 4)

  expected_aug = tibble(
    final_idx = c(6:11, 7:12, 8:13),
    final_datetime = c(as.Date("2017-01-06") + 0:5, as.Date("2017-01-07") + 0:5, as.Date("2017-01-08") + 0:5),
    shift = c(rep(-1, 6), rep(0, 6), rep(1, 6))
  )
  expected_fragments = bind_cols(bind_rows(d, d, d), expected_aug)
  expected_losses = tibble(
    k = c(1,1,1),
    shift = c(-1,0,1),
    loss = c(3,4,5)
  )

  out = data %>% anneal(
    fitted_obs = fitted_obs,
    digest = d,
    range_start = -1,
    range_end = 1,
    loss_fn = rmse
    )
  expect_equal(out$losses, expected_losses)
  expect_equal(out$fragments, expected_fragments)
})


test_that("anneal warns no future vals in fragment", {
  data = tibble(
    x = 1:8,
    y = 1:8,
    datetime = as.Date("2017-01-01") + 0:7,
  ) %>%
    as_tsibble()

  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, newdata = data %>% select(x))
  data = data %>% mutate(fitted_obs = pred_lm)
  d = data %>% digest(datetime, y, 2, 4, n_future_steps = 1)

  # Need to nest multi-warning cases.
  expect_warning(
    expect_warning(
      data %>% anneal(
        fitted_obs = fitted_obs,
        digest = d,
        range_start = -2,
        range_end = 1,
        loss_fn = rmse
      )
    )
  )
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
    adj_idx = c(11:16, 11:20),
    adj_datetime = c(as.Date("2017-01-11") + 0:5, as.Date("2017-01-11") + 0:9)
  )
  out = data %>% digest(datetime, x, 2, 4)
  expect_equal(out, expected)
})

test_that("digest, no partial overlaps", {
  data = tibble(
    x = 101:112,
    datetime = as.Date("2017-01-01") + 0:11
  )
  expected = tibble(
    x = c(107:112, 103:112),
    datetime = c(as.Date("2017-01-07") + 0:5, as.Date("2017-01-03") + 0:9),
    idx = c(7:12, 3:12),
    k = c(rep(1,6), rep(2,10)),
    adj_idx = c(11:16, 11:20),
    adj_datetime = c(as.Date("2017-01-11") + 0:5, as.Date("2017-01-11") + 0:9)
  )
  out = data %>% digest(datetime, x, 2, 4)
  expect_equal(out, expected)
})

test_that("digest, partial overlaps included", {
  data = tibble(
    x = 101:112,
    datetime = as.Date("2017-01-01") + 0:11
  )
  expected = tibble(
    x = c(104:112, 101:112),
    datetime = c(as.Date("2017-01-04") + 0:8, as.Date("2017-01-01") + 0:11),
    idx = c(4:12, 1:12),
    k = c(rep(1,9), rep(2,12)),
    adj_idx = c(8:16, 9:20),
    adj_datetime = c(as.Date("2017-01-08") + 0:8, as.Date("2017-01-09") + 0:11)
  )
  out = data %>% digest(datetime, x, 5, 4)
  expect_equal(out, expected)
})

test_that("digest, n_future_steps limits fragment size", {
  data = tibble(
    x = 101:112,
    datetime = as.Date("2017-01-01") + 0:11
  )
  expected = tibble(
    x = c(104:111, 101:107),
    datetime = c(as.Date("2017-01-04") + 0:7, as.Date("2017-01-01") + 0:6),
    idx = c(4:11, 1:7),
    k = c(rep(1,8), rep(2,7)),
    adj_idx = c(8:15, 9:15),
    adj_datetime = c(as.Date("2017-01-08") + 0:7, as.Date("2017-01-09") + 0:6),
  )
  out = data %>% digest(datetime, x, 5, 4, n_future_steps = 3)
  expect_equal(out, expected)
})

test_that("digest, max_na_sequence", {
  data = tibble(
    x = c(101,102,NA,104,NA,NA,107,NA,NA,NA,NA,112),
    datetime = as.Date("2017-01-01") + 0:11
  )
  expected = tibble(
    x = c(c(107,112),c(102,NA,104,107,112)),
    datetime = c(
      c(as.Date("2017-01-07"), as.Date("2017-01-12")),
      c(as.Date("2017-01-02") + 0:2, as.Date("2017-01-07"), as.Date("2017-01-12"))
    ),
    idx = c(c(7,12),c(2,3,4,7,12)),
    k = c(rep(1,2), rep(2,5)),
    adj_idx = c(c(11,16),c(10,11,12,15,20)),
    adj_datetime = c(
      c(as.Date("2017-01-11"), as.Date("2017-01-16")),
      c(as.Date("2017-01-10") + 0:2, as.Date("2017-01-15"), as.Date("2017-01-20"))
    )
  )
  # k = 1 spans idx 6-12
  # k = 2 spans idx 2-12
  out = data %>% digest(datetime, x, 3, 4, max_na_sequence = 1)
  expect_equal(out, expected)
})


test_that("digest, min_frag_len", {
  data = tibble(
    y = 1:9,
    x = 1:9,
    datetime = as.Date("2017-01-01") + 0:8,
  ) %>%
    as_tsibble()
  expected = tibble(
    y = 4:6,
    x = 4:6,
    datetime = as.Date("2017-01-04") + 0:2,
    idx = 4:6,
    k = 1,
    adj_idx = 8:10,
    adj_datetime = as.Date("2017-01-08") + 0:2
  )

  out = data %>% digest(
    datetime, y, 2, 4,
    n_future_steps = 1,
    min_fragment_len = 3)
  out2 = data %>% digest(
    datetime, y, 2, 4,
    n_future_steps = 1,
    min_fragment_len = 4)

  expect_equal(out, expected)
  expect_equal(nrow(out2), 0)
})

# an internal method (currently)
test_that("mark_long_na_sequences", {
  data = tibble(x = c(NA,3,3,3,NA,NA,2,2,NA,NA,NA,1,NA,2,2,NA,NA))
  expected = c(rep(FALSE,4),rep(TRUE,2),
               rep(FALSE,2),rep(TRUE,3),
               rep(FALSE,4),rep(TRUE,2))
  out = mark_long_na_sequences(data$x, 1)
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

test_that("fragment trimming removes empty fragments", {
  data = tibble(
    x = c(NA,NA,2,2,NA,NA,NA,NA,1,4,NA),
    k = c(1,1,1,1,1,2,2,3,3,3,3)
  )
  expected = tibble(
    x = c(2,2,1,4),
    k = c(1,1,3,3)
  )
  expected_neither = tibble(
    x = c(NA,NA,2,2,NA,NA,1,4,NA),
    k = c(1,1,1,1,1,3,3,3,3)
  )

  expect_equal(data %>% trim_fragments_na(x), expected)
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
