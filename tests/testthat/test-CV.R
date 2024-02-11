test_that("cv anneal output contents", {
  # data fitting setup
  data = tibble(
    x = 1:13,
    y = 1:13,
    datetime = as.Date("2017-01-01") + 0:12,
  ) %>%
    as_tsibble()
  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, new_data = data$x)
  data = data %>% mutate(fitted_obs = pred_lm)

  # args setup
  outer_digest_args = list(n_overlap = 4, n_future_steps = 3)
  inner_digest_args = list(n_overlap = 3, n_future_steps = 2)
  anneal_args = list(range_start = 0, range_end = 0)
  season_len = 4

  # fn call
  cv_out = data %>% anneal_cv_anneal(
    datetime,
    y,
    fitted_obs,
    season_len = season_len,
    outer_digest_args = outer_digest_args,
    inner_digest_args = inner_digest_args,
    anneal_args = anneal_args
  )

  # args come out backend
  expect_equal("fitted_obs", cv_out$.smoothed_obs_name)
  expect_equal(outer_digest_args$n_overlap, cv_out$outer_digest_args$n_overlap)
  expect_equal(outer_digest_args$n_future_steps, cv_out$outer_digest_args$n_future_steps)
  expect_equal(inner_digest_args$n_overlap, cv_out$inner_digest_args$n_overlap)
  expect_equal(inner_digest_args$n_future_steps, cv_out$inner_digest_args$n_future_steps)
  expect_equal(anneal_args$range_start, cv_out$anneal_args$range_start)
  expect_equal(anneal_args$range_end, cv_out$anneal_args$range_end)
  expect_equal(season_len, cv_out$outer_digest_args$season_len)
  expect_equal(season_len, cv_out$inner_digest_args$season_len)

  # should get four-fold cv (four seasons)
  expect_equal(4, length(cv_out$anneals))
})

test_that("cv anneal annealing", {
  # data fitting setup
  data = tibble(
    x = 1:13,
    y = 1:13,
    datetime = as.Date("2017-01-01") + 0:12,
  ) %>%
    as_tsibble()
  fit_lm = lm(y ~ x, data = data)
  pred_lm = predict(fit_lm, new_data = data$x)
  data = data %>% mutate(fitted_obs = pred_lm)

  # args setup
  outer_digest_args = list(n_overlap = 4, n_future_steps = 3)
  inner_digest_args = list(n_overlap = 3, n_future_steps = 2)
  anneal_args = list(range_start = 0, range_end = 0)
  season_len = 4

  # fn call
  cv_out = data %>% anneal_cv_anneal(
    datetime,
    y,
    fitted_obs,
    season_len = season_len,
    outer_digest_args = outer_digest_args,
    inner_digest_args = inner_digest_args,
    anneal_args = anneal_args
  )

  # Looking at the losses should get you at least mostly there.
  losses_1 = tibble(
    k = c(2,3,4),
    shift = 0,
    loss = c(4,8,NaN)
  )
  losses_2 = tibble(
    k = c(1,3,4),
    shift = 0,
    loss = c(4,4,8)
  )
  losses_3 = tibble(
    k = c(1,2,4),
    shift = 0,
    loss = c(8,4,4)
  )
  losses_4 = tibble(
    k = c(1,2,3),
    shift = 0,
    loss = c(12,8,4)
  )
  expect_equal(losses_1, cv_out[["anneals"]][["1"]]$anneal_output$losses)
  expect_equal(losses_2, cv_out[["anneals"]][["2"]]$anneal_output$losses)
  expect_equal(losses_3, cv_out[["anneals"]][["3"]]$anneal_output$losses)
  expect_equal(losses_4, cv_out[["anneals"]][["4"]]$anneal_output$losses)
})
