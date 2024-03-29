#' Anneal fragments of a digest to your data.
#'
#' @param data A tsibble. Cannot have multiple elements in key, if key exists.
#' @param fitted_obs The column in `data` (assumed same name in `digest`)
#'   of fitted observations from `predict(loess_fit, new_data = <your data idx>)`.
#' @param digest Output of `anneal::digest` for smoothed data.
#' @param range_start Start of range to search over.
#' @param range_end End of range to search over.
#' @param loss_fn f(original, fragment) The loss function.
#' @returns list(
#'   shifted_data = tibble. Your original data, repeated (range_end - range_start)
#'                  times, with each repetition augmented by:
#'     final_idx       = The shifted index as aligned to the original data's indexes.
#'     final_datetime  = The `final_idx` in terms of the original data's datetime index.
#'     shift           = The shift along x-axis, relative you the `digest`'s `adj_idx`.
#'                       `adj_idx` + `shift` = `final_idx`
#'   ,
#'   losses       = tibble[
#'     k         = The fragment (1 = 1 season back, etc.).
#'     shift     = The shift along x-axis, relative you the `digest`'s `adj_idx`.
#'                 `adj_idx` + `shift` = `final_idx`
#'     loss      = The loss in the overlapping region of the fragment at `final_idx`
#'                 and the `data`, as calculated by `loss_fn`.
#'   ])
#'
#' @export
anneal <- function(data, fitted_obs, digest, range_start, range_end, loss_fn) {
  data = data %>% mutate(idx_data = row_number())
  # TODO: drop loss col from digest

  df = NULL
  losses = NULL
  for (k_idx in (digest %>% pull(k) %>% unique())) {
    fragment = digest %>% filter(k == k_idx)
    annealed_frag = data %>% anneal_fragment({{ fitted_obs }}, fragment, range_start, range_end, loss_fn)
    df = bind_rows(df, annealed_frag$fragments)
    losses = bind_rows(losses, annealed_frag$losses)
  }

  list(fragments = df, losses = losses)
}


anneal_fragment = function(data, fitted_obs, fragment, range_start, range_end, loss_fn) {
  df = NULL
  losses = NULL
  for (shift in seq(range_start, range_end)) {
    shifted_fragment = fragment %>%
      mutate(
        final_idx = adj_idx + shift,
        final_datetime = adj_datetime + shift,
        shift = shift
      )

    # how far beyond the "real" ts does our shifted fragment extend?
    # n_to_append = max(
    #   0,
    #   max(shifted_fragment$final_idx) - max(data$idx_data)
    # )  # If <= 0 you can't forecast.
    # if (n_to_append == 0) {
    #   frag_k = fragment %>% pull(k) %>% first()
    # }

    # Loss calc
    overlap = inner_join(
      data,
      shifted_fragment,
      by = join_by(idx_data == final_idx)
    )  # -> tsibble
    # TODO: check overlap [min,max]
    fitted_obs_col_name = data %>%
      as_tibble() %>%
      ungroup() %>%
      select({{ fitted_obs }}) %>%
      names()
    fitted_obs_x_name = paste0(fitted_obs_col_name, ".x")
    fitted_obs_y_name = paste0(fitted_obs_col_name, ".y")
    loss = loss_fn(
      overlap %>% pull(fitted_obs_x_name),
      overlap %>% pull(fitted_obs_y_name)
    )
    loss_row = tibble(
      k = fragment %>% pull(k) %>% first(),
      shift = shift,
      loss = loss
    )

    # Write
    df = bind_rows(df, shifted_fragment)
    losses = bind_rows(losses, loss_row)
  }

  list(fragments = df, losses = losses)
}


#' Fragment your time series by seasons, with an overlap.
#'
#' A "fragment" is 1 or more prior seasons' of data. From a baseline
#' of k seasons before your latest observation:
#' 1. It includes a "future" portion of `n_future_steps` observations,
#' for forecasting purposes.
#' 2. It also includes a "past" portion of `n_overlap` observations,
#' for annealing purposes.
#' This function is of the same ilk as `stretch_tsibble` from the `tsibble` package.
#'
#' The data cannot contain gaps, because this function relies on indexing.
#' You can be sure that it doesn't by making it a `tsibble` and checking
#' for gaps.
#'
#' Your `n_future_steps` should be longer than your desired forecast horizon,
#' because `anneal` will potentially shift it backwards. A sensible lower bound
#' is (your desired forecast horizon) + (your max backwards shift in `anneal`).
#'
#' `max_na_sequence` is used to remove long sequences of missing data.
#' Because `predict` can predict anywhere (i.e. the observation can be NA),
#' this prevents them from being `predict`ed in `anneal` and thereby
#' prevents them from factoring into `anneal`'s loss calculation.
#'
#' Fragments shorter than `min_fragment_len`, after any trimming via
#' `max_na_sequence` (though the fragment may still contain or
#' have at its ends or even be entirely NA observations,
#' which count towards the length) are removed.
#' It defaults to 2, because:
#'   - What can you do with a fragment of length 1?
#'   - tsibble doesn't like it.
#'
#' @param data tsibble. The data. Must not contain gaps.
#' @param .datetime Your datetime column. The tsibble's "index".
#' @param .observation Your observations column.
#' @param n_overlap The number of time points to overlap. Must be positive integer.
#' @param season_len The number of time points per season. Must be positive integer.
#' @param n_future_steps The number of time points for fragments to extend beyond your latest observation.
#' @param max_na_sequence The largest series of NAs (i.e. gap between observations) before filtering out.
#' @param min_fragment_len The minimum length of a fragment to include in output.
#' @returns A tibble of fragments augmented with cols:
#'   idx:          The index w.r.t the original data.
#'   k:            The fragment (1 = 1 season back, etc.).
#'   adj_idx:      The fragment's index aligned to the latest data.
#'   adj_datetime: The fragment's datetime aligned to the latest data.
#'
#' @export
digest <- function(data,
                   .datetime,
                   .observation,
                   n_overlap,
                   season_len,
                   n_future_steps = 120,
                   max_na_sequence = Inf,
                   min_fragment_len = 2) {
  assert_that(n_overlap > 0, msg = "Need n_overlap > 0.")
  assert_that(season_len > 0, msg = "Need season_len > 0.")
  assert_that(n_future_steps > 0, msg = "Need n_future_steps > 0.")

  data = data %>%
    as_tibble() %>%
    ungroup() %>%
    mutate(idx = row_number())

  data = data %>%
    mutate(to_remove = mark_long_na_sequences({{ .observation }}, max_na_sequence))

  # Remember: R indexes start at 1.
  k = 1
  df = NULL
  t_max = data %>% select(idx) %>% max()
  while (TRUE) {
    i = t_max - (k * season_len) + 1  # `+ 1` needed
    if (i <= 1) {
      break  # Nothing to overlap with.
    }

    j = max(i - n_overlap, 1)
    fragment = data %>%
      slice(j:n()) %>%
      mutate(
        k = k,
        adj_idx = idx + (k * season_len),
        adj_datetime = {{ .datetime }} + (k * season_len)
      )
    df = bind_rows(df, fragment)
    k = k + 1
  }
  df = df %>% filter(adj_idx <= t_max + n_future_steps)
  # Removing large sequences of NA
  df = df %>%
    filter(!to_remove) %>%  # note the `!`
    select(-to_remove)
  # Removing short fragments
  df = df %>%
    group_by(k) %>%
    filter(n() >= min_fragment_len) %>%
    ungroup()

  df
}


#' Create a boolean vector indicating long sequences of NAs.
#'
#' This function is useful for subsequent use in filtering out these long
#' sequences from your data.
#'
#' Note that the long sequences of NAs are the "TRUE"s, so to filter,
#' be mindful of using the logical not `!`.
#'
#' @param observations The observations to generate the boolean vector for.
#' @param maxgap The longest sequence of NAs to lump in (boolean-ly) with the non-NAs.
#'
#' @returns A boolean vector of length(observations).
#'   TRUE  = The observation (well, lack thereof) is part of a sequence of NAs
#'           that exceeds maxgap in length.
#'   FALSE = Either observations, or sequences of NAs <= maxgap.
mark_long_na_sequences <- function(observations, maxgap) {
  # ref: https://github.com/SteffenMoritz/imputeTS/blob/a3155b041f1a7d36bdec3c152043a90c948f225e/R/na_seasplit.R#L210

  # Get logical vector of the time series via is.na() and then get the
  # run-length encoding of it. The run-length encoding describes how long
  # the runs of FALSE (is not NA) and TRUE (is NA) are
  rlencoding <- rle(is.na(observations))

  # Runs smaller than or equal to maxgap are set FALSE
  rlencoding$values[rlencoding$lengths <= maxgap] <- FALSE

  # The original vector is being reconstructed by inverse.rls, only now the
  # shorter runs of NAs are, like the non-NA observations, assigned "FALSE".
  # Only the runs of NAs exceeding `maxgap` are "TRUE".
  inverse.rle(rlencoding)
}

#' Trim digest fragments of leading/trailing NAs
#'
#' This is useful for cases where one/both ends of your fragment(s)'s
#' original observations contain NA's that you don't want to have `predict`ed with
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
