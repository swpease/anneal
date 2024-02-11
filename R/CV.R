#' Generate (# of seasons)-fold CV data.
#'
#' Generate CV data using a particular set of hyperparameters.
#'
#' This function may be plugged into a grid search for hyperparameters.
#' The `.datetime` and `.observation` columns should stay the
#' same across calls. Most of the other args could reasonably be CV'ed. A
#' few things worth mentioning:
#'
#' `tmax_offset` is useful for seeing how a specific method / set of hyperparams
#' performs at different times of the season. For instance, for daily data, you
#' could use values such as -14, -7, 0, 7, 14, to help evaluate +- two weeks.
#'
#' The `loss_fn` arg for `anneal`, which you would pass via `anneal_args`,
#' needs to be a list (e.g. list(myfun = myfun)) for reasons.
#'
#' See the source code for the default args in the o/i/a-args lists.
#'
#' The rationale behind this function is that the smoothed data was locally
#' smoothed, so CV should just work okay on a seasonal basis -- Jan 13 of
#' this year has no bearing on Jan 13 of next year, w.r.t smoothing. Ergo,
#' there is a two-step digest process:
#'   1. An outer, broad-range digest, yielding a set of k (#-seasons of) "data"s.
#'   2. An inner, more restricted digest, yielding a set of k (#-seasons of)
#'      "fragments".
#' On each interation i, the outer "data" fragment k == i is used, and the
#' inner fragments k != i are used.
#' Each outer fragment is padded back to t_0 with NAs for compatibility with the
#' current implementation of `anneal`.
#'
#' I've padded the data with another season of NAs. That means that the k == 1
#' fragment is the current season, while in a typical `digest`, k == 1 is prior
#' season (because the current season isn't even included). This should yield
#' a bit more power when `tmax_offset` < 0 and/or the annealed k == 1 fragment
#' has a positive offset, and I couldn't think of a downside. I may make it
#' optional, though.
#'
#' @param data tsibble. Your data, with a column of smoothed observations.
#' @param .datetime `data`'s column of datetimes, i.e. your tsibble's `index`.
#' @param .observation `data`'s column of original observations.
#' @param .smoothed_obs `data`'s column of smoothed observations.
#' @param tmax_offset int. The offset from your latest observation to consider
#'                    as "tmax".
#' @param season_len int. The season length.
#' @param outer_digest_args list. Any arguments to override the defaults for
#'                          the outer digest loop.
#' @param inner_digest_args list. Any arguments to override the defaults for
#'                          the inner digest loop.
#' @param anneal_args list. Any arguments to override the defaults for `anneal`.
#' @returns list. All of the arguments used from annealing, i.e.:
#'            - .smoothed_obs_name,
#'            - tmax_offset,
#'            - (outer/inner/anneal)-args
#'          Also, annealing information. A list, with an element for each outer-
#'          loop digest fragment.
#'            - anneal_output,
#'            - test_data_outer_fragment (i.e. the outer-loop digest fragment),
#'            - k_test, i.e. the k of the test_data_outer_fragment
#'
#' @export
anneal_cv_anneal <- function(data,
                             .datetime,
                             .observation,
                             .smoothed_obs,
                             tmax_offset = 0,
                             season_len = 365,
                             outer_digest_args = list(),
                             inner_digest_args = list(),
                             anneal_args = list()) {
  # Digest args setup
  default_outer_digest_args = list(
    n_overlap = 180,
    season_len = season_len,  # == 365
    n_future_steps = 120,
    max_na_sequence = Inf  # Bounds here ought to be redundant w/ inner digest.
    # ...b/c loss fn will have NA's on inner digest's fragments, you'll just
    # get those indexes omitted via that anyway.
  )
  default_inner_digest_args = list(
    n_overlap = 60,
    season_len = season_len,  # == 365
    n_future_steps = 80,
    max_na_sequence = 100
  )
  od_args = modifyList(default_outer_digest_args, outer_digest_args)
  id_args = modifyList(default_inner_digest_args, inner_digest_args)

  # Anneal args setup
  default_anneal_args = list(
    range_start = -20,
    range_end = 20,
    loss_fn = list(slope_weighted_rmse = slope_weighted_rmse)
  )
  anneal_args = modifyList(default_anneal_args, anneal_args)

  # Getting .smoothed_obs col name
  smoothed_obs_col_name = data %>%
    as_tibble() %>%
    ungroup() %>%
    select({{ .smoothed_obs }}) %>%
    names()

  # Add rows up to a year from your desired "now" (i.e. your most recent
  # observation, adjusted by any offset per `tmax_offset`).
  # This lets the current season be used (to at least some extent)
  # in CV. It only really helps for a negative `tmax_offset`.
  # TODO: do I want to make this an option?
  data = data %>% append_row(n = (season_len + tmax_offset))

  outer_digest = data %>%
    digest(
      {{ .datetime }},
      {{ .observation }},
      n_overlap = od_args$n_overlap,
      season_len = od_args$season_len,
      n_future_steps = od_args$n_future_steps,
      max_na_sequence = od_args$max_na_sequence
    )
  inner_digest = data %>%
    digest(
      {{ .datetime }},
      {{ .observation }},
      n_overlap = id_args$n_overlap,
      season_len = id_args$season_len,
      n_future_steps = id_args$n_future_steps,
      max_na_sequence = id_args$max_na_sequence
    )

  # Set-up list for the k CV outputs
  ks = outer_digest %>% pull(k) %>% unique()
  anneals = vector(mode = "list", length = length(ks))
  names(anneals) <- ks
  # The entire output object. Woof!
  output_info = list(
    .smoothed_obs_name = smoothed_obs_col_name,
    tmax_offset = tmax_offset,
    outer_digest_args = od_args,
    inner_digest_args = id_args,
    anneal_args = anneal_args,
    anneals = anneals
  )

  # The actual CV annealing.
  for (k_idx in ks) {
    # Anneal expect a tsibble.
    # TODO: Is it's tsibble-ness used ever anymore?
    # In any case, need to pad back to t=1 from the adj_idx,
    # b/c `anneal` sets-up its own index using `row_number`.
    test_frag = outer_digest %>%
      filter(k == k_idx) %>%
      as_tsibble(index = adj_idx)
    min_idx = test_frag %>% select(adj_idx) %>% min()
    test_frag = test_frag %>% append_row(-(min_idx - 1))

    train_frags = inner_digest %>% filter(k != k_idx)
    out = test_frag %>%
      anneal(
        {{ .smoothed_obs }},
        train_frags,
        anneal_args$range_start,
        anneal_args$range_end,
        anneal_args$loss_fn[[1]]
      )
    output_info$anneals[[k_idx]] = list(
      anneal_output = out,
      test_data_outer_fragment = test_frag,
      k_test = k_idx
    )
  }

  output_info
}
