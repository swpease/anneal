#' Plot the output of anneal. Tentative.
#'
#' @export
plot_anneal <- function(data, x, y, anneal_output) {
  # anneal_output$losses %>%
  #   ggplot() +
  #   geom_point(aes(shift, loss))

  ggplot() +
    geom_line(
      data = anneal_output$fragments,
      mapping = aes(
        x = final_idx,
        y = .pred_obs,
        group = shift
      ),
      color = "blue"
    ) +
    geom_line(
      data = data,
      mapping = aes(
        x = {{ x }},
        y = {{ y }}
      )
    )
}


#'
#' Assumptions:
#'   1. Object passed to `loess_fit` had the formula `[var1] ~ [var2]`, e.g. `obs ~ date_time`.
#'   2. Object passed to `data` includes the `[var1]` col from Assumption #1 with the same name, e.g. `obs`.
#'
#' @param data A tsibble. Cannot have multiple elements in key, if key exists.
#' @param fitted_obs_col_name The column name in `data` of fitted observations from `predict(loess_fit, new_data = <your data idx>)`.
#' @param digest Output of `anneal::digest` for smoothed data.
#' @param resolution The step size (along x-axis) during search. Currently must equally divide 1 (e.g. 1, 0.2, or 0.5).
#' @param range_start Start of range to search over.
#' @param range_end End of range to search over.
#' @param loess_fit Fitted output of `loess` on your data.
#' @param loss_fn f(original, fragment) The loss function.
#' @returns list(
#'   shifted_data = tibble[
#'     k         = The fragment (1 = 1 season back, etc.).
#'     .pred_obs = The fragment's predicted value at `final_idx` for a particular `shift`.
#'     final_idx = The index as aligned to the original data's indexes.
#'     shift     = The shift along x-axis, relative you the `digest`'s `adj_idx`.
#'                 `adj_idx` + `shift` = `final_idx`
#'     datetime  = The `final_idx` in terms of the original data's datetime index.
#'   ],
#'   losses       = tibble[
#'     shift     = The shift along x-axis, relative you the `digest`'s `adj_idx`.
#'                 `adj_idx` + `shift` = `final_idx`
#'     loss      = The loss in the overlapping region of the fragment at `final_idx`
#'                 and the `data`, as calculated by `loss_fn`.
#'   ])
#'
#' @export
anneal <- function(data, fitted_obs_col_name, digest, resolution, range_start, range_end, loess_fit, loss_fn) {
  data = data %>% mutate(idx_data = row_number())
  # TODO: drop loss col from digest

  df = NULL
  losses = NULL
  for (k_idx in (digest %>% pull(k) %>% unique())) {
    fragment = digest %>% filter(k == k_idx)
    annealed_frag = anneal_fragment(data, fitted_obs_col_name, fragment, resolution, range_start, range_end, loess_fit, loss_fn)
    df = bind_rows(df, annealed_frag$fragments)
    losses = bind_rows(losses, annealed_frag$losses)
  }

  list(fragments = df, losses = losses)
}


anneal_fragment = function(data, fitted_obs_col_name, fragment, resolution, range_start, range_end, loess_fit, loss_fn) {
  pred_col_name = attr(loess_fit$terms, "term.labels")
  upsampled_fragment = tibble(
    idx_upsam = seq(min(fragment$idx), max(fragment$idx), resolution),
    adj_idx = seq(min(fragment$adj_idx), max(fragment$adj_idx), resolution),
    k = fragment %>% pull(k) %>% first(),
    .pred_obs = predict(loess_fit, tibble({{ pred_col_name }} := idx_upsam))
  )

  df = NULL
  losses = NULL
  for (shift in seq(range_start, range_end, resolution)) {
    # shift
    shifted_upsampled_fragment = upsampled_fragment %>%
      mutate(
        final_idx = adj_idx + shift,
        shift = shift
      )

    # Remove non-(near)-integer indexes
    shifted_upsampled_fragment <- shifted_upsampled_fragment %>%
      filter(near(final_idx %% 1, 1) | near(final_idx %% 1, 0)) %>%  # `near` in case of, e.g. res = 0.3333
      mutate(final_idx = round(final_idx))

    # Add date col.
    # how far beyond the "real" ts does our shifted fragment extend?
    n_to_append = max(shifted_upsampled_fragment$final_idx) - max(data$idx_data)
    # TODO: include this check?
    # assert_that(n_to_append > 0, msg = "Fragment doesn't extend beyond original series.")
    # extend our dates that far
    extended_dates = data %>%
      select(index(.)) %>%  # ref: https://magrittr.tidyverse.org/reference/pipe.html#using-the-dot-for-secondary-purposes
      append_row(n = n_to_append) %>%
      mutate(idx_dates = row_number())  # -> tsibble[date_col]
    # join
    shifted_upsampled_fragment <- left_join(
      shifted_upsampled_fragment,
      extended_dates,
      by = join_by(final_idx == idx_dates),
      relationship = "one-to-one"
    )  # date col will be whatever the `index` of the input `data` tsibble is.

    # Loss calc
    overlap = inner_join(
      data,
      shifted_upsampled_fragment,
      by = join_by(idx_data == final_idx)
    )  # -> tsibble
    # TODO: check overlap [min,max]
    original_obs = overlap %>%
      as_tibble() %>%
      ungroup() %>%
      pull(fitted_obs_col_name)  # woof
    loss = loss_fn(original_obs, overlap$.pred_obs)
    loss_row = tibble(
      k = fragment %>% pull(k) %>% first(),
      shift = shift,
      loss = loss
    )

    # Write
    df = bind_rows(df, shifted_upsampled_fragment)
    losses = bind_rows(losses, loss_row)
  }

  list(fragments = df, losses = losses)
}


#' Calculate the root mean squared error of two sequences.
#'
#' @param a Sequence 1
#' @param b Sequence 2
#' @returns RMSE
#'
#' @export
rmse <- function(a, b) {
  sqrt(mean((a - b) ^ 2))
}


#' Calculate the weighted root mean squared error of two sequences.
#'
#' Weights are based on sequence `b`'s `difference`. At time t,
#' the formulation is: ((a_t - b_t) * (b_t - b_{t-1})) ^ 2.
#'
#' Intended to be passed as the `loss_fn` arg of `anneal`.
#'
#' In the context of this package, `b` makes more sense as your fragment,
#' because the "current" data, `a`, may have a misleading trajectory for
#' its smoothed fit (e.g. an imminent upturn / other inflection point
#' isn't reflected because it hasn't happened yet.).
#'
#' This weighting means that the flat parts / inflection points of
#' the series have minimal influence. Ergo, the sequences' alignment
#' should be driven by the spring/fall temperature changes, not the
#' extremeness of summer/winter. It is a compromise for the failure
#' of the method relying on calculating an orthogonal vector.
#'
#' @param a Sequence 1
#' @param b Sequence 2
#' @returns Weighted RMSE
#'
#' @export
weighted_rmse <- function(a, b) {
  b_diff = difference(b)

  sqrt(mean(((a - b) * b_diff) ^ 2, na.rm = TRUE))
}


#' Fragment your time series by seasons, with an overlap.
#'
#' A "fragment" is 1 or more prior seasons' of data, plus an overlap for annealing.
#' This function is of the same ilk as `stretch_tsibble` from the `tsibble` package.
#'
#' Note that potentially "empty" fragments may be returned (e.g. no data
#' collected for a given year). To remove these, you can pass the output to
#' `trim_fragments_na`.
#'
#' @param data The data.
#' @param n_overlap The number of time points to overlap. Must be positive integer.
#' @param season_len The number of time points per season. Must be positive integer.
#' @param n_future_steps The number of time points for fragments to extend beyond your latest observation.
#' @param include_partial_overlap Whether to include the final fragment if there is only partial overlap.
#' @returns A tibble of fragments containing cols:
#'   idx:     The index w.r.t the original data.
#'   k:       The fragment (1 = 1 season back, etc.).
#'   adj_idx: The fragment, aligned to the latest data.
#'
#' @export
digest <- function(data, n_overlap, season_len, n_future_steps = 120, include_partial_overlap = TRUE) {
  assert_that(n_overlap > 0, msg = "Need n_overlap > 0.")
  assert_that(season_len > 0, msg = "Need season_len > 0.")
  assert_that(n_future_steps > 0, msg = "Need n_future_steps > 0.")

  data = data %>%
    as_tibble() %>%
    ungroup() %>%
    mutate(idx = row_number())

  # Remember: R indexes start at 1.
  k = 1
  df = NULL
  t_max = data %>% select(idx) %>% max()
  while (TRUE) {
    i = t_max - (k * season_len) + 1  # `+ 1` needed
    if (i <= 1) {
      break  # Nothing to overlap with.
    }

    j = i - n_overlap
    if (j <= 0) {
      if (!include_partial_overlap) {
        break
      } else {
        sprintf("Incomplete overlap for longest fragment. Missing %s of %s values.", (-j + 1), n_overlap)
      }
    }

    j = max(j, 1)
    fragment = data %>%
      slice(j:n()) %>%
      mutate(
        k = k,
        adj_idx = idx + (k * season_len)
      )
    df = bind_rows(df, fragment)
    k = k + 1
  }
  df = df %>% filter(adj_idx <= t_max + n_future_steps)

  df
}


#' Trim digest fragments of leading/trailing NAs
#'
#' This is useful for cases where one/both ends of your fragments contain
#' NA's that you don't want to have `predict`ed with
#' your `loess_fit` in `anneal`, such as a large gap beyond the end(s) of the
#' fragment until the next non-NA observation (i.e. you don't trust the fit
#' for these NA's).
#'
#' This function always removes fragments containing only NA observations.
#' As such, setting both `left` and `right` to `FALSE` will remove
#' only these all-NA fragments.
#'
#' @param digest Output of `digest`.
#' @param .col The column to search through for NAs.
#' @param left Trim left side?
#' @param right Trim right side?
#' @returns digest with trimmed fragments.
#'
#' @export
trim_fragments_na <- function(digest, .col, left = TRUE, right = TRUE) {
  # Remove fragments w/o any non-NA observations.
  digest = digest %>%
    group_by(k) %>%
    filter(!all(is.na({{ .col }}))) %>%
    ungroup()

  if (left) {
    digest = digest %>% slice(
      detect_index(
        {{ .col }},
        \(x) !is.na(x)
      ):n(),
      .by = k)
  }
  if (right) {
    digest = digest %>% slice(
      1:detect_index(
        {{ .col }},
        \(x) !is.na(x),
        .dir = "backward"
      ),
      .by = k)
  }

  digest
}


#' Calculate a vector orthogonal to the normalized average difference of your time series.
#'
#' The steps are:
#'   1. Take the last n values of your ts.
#'   2. Take their difference.
#'   3. Perform a 90 degree counter-clockwise rotation on (diff, 1) vectors.
#'   4. Normalize each rotated vector.
#'   5. Optionally weight each vector.
#'   6. Add these vectors.
#'   7. Normalize the total vector.
#'
#' @param ts The time series. Should be smooth/smoothed. e.g. dat$readings
#' @param n The number of values at the end of the value column to include. Constraint: 1 < n < len(value)
#' @param scaler Scale your y's so that they're on equal footing w/ dx.
#' @param wts Optional weights for each difference.
#'
#' @returns tibble (1x2) Normalized orthogonal vector.
#'
#' @export
calc_ortho_vec <- function(ts, n, scaler = 1, wts = NULL) {
  assert_that(n > 1, msg = "Need n > 1.")
  validate_that(n < length(ts), msg = "n >= len(value); will yield NAs for diff calc.")
  if (!(is.null(wts))) {
    validate_that(length(n) == length(wts), msg = "Want one wt per difference.")
  }
  wts = ifelse(is.null(wts), 1, wts)

  diffs = ts %>% difference() %>% as_tibble() %>% slice_tail(n = n)
  scaled_diffs = diffs * scaler
  rotated_vecs = tibble(x = -scaled_diffs, y = 1)  # 90deg counter-clock rot.
  normed_rotated_vecs = rotated_vecs %>%  # TODO: extract
    mutate(
      vec_size = sqrt(x^2 + y^2),
      x = x / vec_size,
      y = y / vec_size
    ) %>%
    select(-vec_size)
  weighted_normed_rotated_vecs = normed_rotated_vecs %>%
    mutate(
      x = x * wts,
      y = y * wts
    )
  total_vec = weighted_normed_rotated_vecs %>%
    summarise(
      x = sum(x),
      y = sum(y)
    )
  normed_total_vec = total_vec %>%
    mutate(
      vec_size = sqrt(x^2 + y^2),
      x = x / vec_size,
      y = y / vec_size
    ) %>%
    select(-vec_size)

  normed_total_vec
}

#' Plot the output of `calc_ortho_vec`.
#'
#' @param normed_total_vec Output of `calc_ortho_vec`.
#'
#' @export
plot_ortho_vec <- function(normed_total_vec) {
  normed_total_vec = normed_total_vec %>% bind_rows(-normed_total_vec)  # Add vec in opposite dir.
  normed_total_vec %>%
    ggplot() +
    geom_segment(
      aes(x = 0, y = 0, xend = x, yend = y, color = "red"),
      arrow = arrow(),
      show.legend = FALSE
    ) +
    xlim(-1.2, 1.2) +
    ylim(-1.2, 1.2)
}

#' Plot a smoothed time series with the orthogonal vector to its head overlaid.
#'
#' The orthogonal vector comes from `calc_ortho_vec`.
#' The smoothed time series is, e.g. a `lowess` output.
#'
#' @param data Your data containing a smoothed time series.
#' @param index The date/time column.
#' @param value The column of smoothed values.
#' @param n The `n` used in `calc_ortho_vec`.
#' @param normed_total_vec Output of `calc_ortho_vec`.
#' @param arrow_scale Scale of arrow. `1` = 1 unit of index, 1 unit of value.
#'
#' @export
plot_data_with_ortho_vec <- function(data, index, value, n, normed_total_vec, arrow_scale = 1) {
  data = data %>% as_tibble() %>% ungroup()  # tsibble keeps idx and keys linked.
  data = data %>% mutate({{ index }} := as_datetime({{ index }}))
  ortho_pts = data %>% slice_tail(n = n)
  mid_pt = data %>% slice((n() - floor(n / 2)))

  # Get "one" (diff b/w time units of index)
  a = data %>% select({{ index }}) %>% slice(n() - 1) %>% pull()
  b = data %>% select({{ index }}) %>% slice(n()) %>% pull()
  dx = as.duration(lubridate::interval(a, b))
  # Starting point of arrow vec
  x0 = mid_pt %>% select({{ index }}) %>% pull()
  y0 = mid_pt %>% select({{ value }}) %>% pull()
  # Point arrow both ways.
  normed_total_vec = normed_total_vec %>% bind_rows(-normed_total_vec)
  # Resize stem len
  normed_total_vec = normed_total_vec * arrow_scale

  ggplot() +
    geom_line(
      data = data,
      aes(x = {{ index }}, y = {{ value }})
    ) +
    geom_line(
      data = ortho_pts,
      aes(x = {{ index }}, y = {{ value }}),
      show.legend = FALSE,
      color = "blue"
    ) +
    geom_segment(
      data = normed_total_vec,
      aes(
        x = x0,
        y = y0,
        xend = (x0 + (x * dx)),
        yend = (y0 + y)
      ),
      arrow = arrow(length = unit(arrow_scale, "points")),
      show.legend = FALSE,
      color = "red"
    )
}
