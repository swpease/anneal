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
  anneal_args = list(
    range_start = 0,
    range_end = 0,
    loss_fn = list("rmse" = rmse)
  )
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
  expect_equal(cv_out$anneal_args$range_start, anneal_args$range_start)
  expect_equal(cv_out$anneal_args$range_end, anneal_args$range_end)
  expect_equal(cv_out$anneal_args$loss_fn, list("rmse" = rmse))
  expect_equal(season_len, cv_out$outer_digest_args$season_len)
  expect_equal(season_len, cv_out$inner_digest_args$season_len)
  expect_equal(cv_out$tmax_idx, 17)


  # should get four-fold cv (four seasons)
  expect_equal(4, length(cv_out$anneals))


  # alt check
  anneal_args2 = list(
    range_start = 0
  )

  # fn call
  cv_out = data %>% anneal_cv_anneal(
    datetime,
    y,
    fitted_obs,
    season_len = season_len,
    outer_digest_args = outer_digest_args,
    inner_digest_args = inner_digest_args,
    anneal_args = anneal_args2
  )
  expect_equal(cv_out$anneal_args$range_start, anneal_args$range_start)
  expect_equal(cv_out$anneal_args$range_end, 20)
  expect_equal(cv_out$anneal_args$loss_fn,
               list(slope_weighted_rmse = slope_weighted_rmse))

})

test_that("cv anneal annealing", {
  # NB: This toy data yields an unusual case of a fragment
  # with only a single observation. tsibble can't handle that,
  # since it doesn't know what the time interval would be,
  # so the current implementation just `next`s these cases.

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
    loss = c(4,4,NaN)
  )
  losses_3 = tibble(
    k = c(1,2,4),
    shift = 0,
    loss = c(8,4,NaN)
  )
  # This fragment is skipped b/c only one element is in the "past",
  # so when it's truncated to not include future values, that one
  # element is all that remains.
  # losses_4 = tibble(
  #   k = c(1,2,3),
  #   shift = 0,
  #   loss = c(12,8,4)
  # )
  expect_equal(losses_1, cv_out[["anneals"]][["1"]]$anneal_output$losses)
  expect_equal(losses_2, cv_out[["anneals"]][["2"]]$anneal_output$losses)
  expect_equal(losses_3, cv_out[["anneals"]][["3"]]$anneal_output$losses)
  # expect_equal(losses_4, cv_out[["anneals"]][["4"]]$anneal_output$losses)
  expect_null(cv_out[["anneals"]][["4"]])
})

test_that("cv anneal warns short fragments", {
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

  expect_warning(data %>% anneal_cv_anneal(
    datetime,
    y,
    fitted_obs,
    season_len = season_len,
    outer_digest_args = outer_digest_args,
    inner_digest_args = inner_digest_args,
    anneal_args = anneal_args
  ))
})
