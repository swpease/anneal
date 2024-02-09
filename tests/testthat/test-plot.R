# ref: https://testthat.r-lib.org/reference/expect_snapshot_file.html#ref-examples

save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  on.exit(dev.off())
  show(code)

  path
}

# ref: https://testthat.r-lib.org/reference/expect_snapshot_file.html#ref-examples
expect_snapshot_plot <- function(name, code) {
  # skips go here
  # -- skips --
  # ref: https://testthat.r-lib.org/articles/skipping.html
  # end skips
  name <- paste0(name, ".png")
  # I don't really get announcing, but it's what they say to do...
  # THIS DOESNT WORK HERE
  announce_snapshot_file(name = name)
  path <- save_png(code)
  expect_snapshot_file(path, name)
}

test_that("end to end", {
  # get data
  sugg_td = readr::read_csv(test_path("SUGG_target_temp_data.csv"))
  td = sugg_td %>%
    as_tsibble()
  # end get data

  # loess
  td = td %>% mutate(idx = row_number())  # loess doesn't like datetimes
  n_loess_pts = 125
  span = n_loess_pts / nrow(td)
  td_loess = loess(
    observation ~ idx,
    td,
    span = span
  )
  smoothed_td = predict(td_loess, newdata = (td %>% select(idx)))
  td = td %>% mutate(smoothed_td = smoothed_td)
  # end loess

  # annealing
  d = td %>% digest(datetime,
                    observation,
                    n_overlap = 55,
                    season_len = 365,
                    max_na_sequence = 100)
  out = td %>%
    anneal(
      fitted_obs = smoothed_td,
      digest = d,
      range_start = -30,
      range_end = 30,
      loss_fn = slope_weighted_rmse
    )
  # end annealing

  # plotting tests
  name = "losses_plot"
  expect_snapshot_plot(name, plot_anneal_losses(out))

  name = "min_loss_fragments_plot"
  expect_snapshot_plot(
    name,
    td %>% plot_anneal_min_loss_fragments(datetime, smoothed_td, out, "smoothed_td")
  )

  name = "original_fragments_plot"
  expect_snapshot_plot(
    name,
    td %>% plot_anneal_original_fragments(datetime, smoothed_td, out, "smoothed_td")
  )

  name = "single_fragment_plot"
  expect_snapshot_plot(
    name,
    td %>% plot_anneal_fragment(datetime, smoothed_td, out, "smoothed_td")
  )
})
