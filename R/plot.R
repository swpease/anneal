#' Plot the min-loss shifted fragments against the original data.
#'
#' Useful for seeing what a forecast would look like,
#' and comparing this to using the naive seasonal via
#' `plot_original_fragments`.
#'
#' Setting t_max at, e.g. `lubridate::today() + 30`, is useful for
#' visual comparisons with `plot_min_loss_fragments` (with the same
#' t_max).
#'
#' @param data The data.
#' @param x The datetime column.
#' @param y The (possibly smoothed) observations.
#' @param anneal_output Output of `anneal`.
#' @param fragment_obs_col_name The (possibly smoothed) fragments' observations column name.
#' @param t_max The upper xlim.
#'
#' @export
plot_anneal_min_loss_fragments <- function(data,
                                           x,
                                           y,
                                           anneal_output,
                                           fragment_obs_col_name,
                                           t_max = NULL) {
  if (is.null(t_max)) {
    t_max = anneal_output$fragments %>%
      pull(final_datetime) %>%
      max()
  }

  min_losses = anneal_output$losses %>%
    group_by(k) %>%
    filter(loss == min(loss, na.rm = TRUE)) %>%
    select(-loss)
  min_loss_fragments = right_join(
    anneal_output$fragments,
    min_losses,
    by = join_by(k, shift)
  )
  min_loss_fragments = min_loss_fragments %>%
    mutate(k = as.factor(k))

  ggplot() +
    geom_line(
      data = min_loss_fragments,
      mapping = aes(
        x = final_datetime,
        y = .data[[fragment_obs_col_name]],
        color = k,
      )
    ) +
    geom_line(
      data = data,
      mapping = aes(
        x = {{ x }},
        y = {{ y }}
      )
    ) +
    scale_x_date(limits = as.Date(c(NA, t_max))) +
    xlab("Datetime") +
    ylab("Predicted Observation") +
    ggtitle("Fragments Shifted to Minimum Loss")
}

#' Plot the unshifted fragments against the original data.
#'
#' Useful for seeing what a naive seasonal forecast would look like,
#' and comparing this to what `anneal` gives you via
#' `plot_min_loss_fragments`.
#'
#' Setting t_max at, e.g. `lubridate::today() + 30`, is useful for
#' visual comparisons with `plot_min_loss_fragments` (with the same
#' t_max).
#'
#' @param data The data.
#' @param x The datetime column.
#' @param y The (possibly smoothed) observations.
#' @param anneal_output Output of `anneal`.
#' @param fragment_obs_col_name The (possibly smoothed) fragments' observations column name.
#' @param t_max The upper xlim.
#'
#' @export
plot_anneal_original_fragments <- function(data,
                                           x,
                                           y,
                                           anneal_output,
                                           fragment_obs_col_name,
                                           t_max = NULL) {
  if (is.null(t_max)) {
    t_max = anneal_output$fragments %>%
      pull(final_datetime) %>%
      max()
  }

  original_fragments = anneal_output$fragments %>%
    filter(shift == 0) %>%
    mutate(k = as.factor(k))

  ggplot() +
    geom_line(
      data = original_fragments,
      mapping = aes(
        x = final_datetime,
        y = .data[[fragment_obs_col_name]],
        color = k,
      )
    ) +
    geom_line(
      data = data,
      mapping = aes(
        x = {{ x }},
        y = {{ y }}
      )
    ) +
    scale_x_date(limits = as.Date(c(NA, t_max))) +
    xlab("Datetime") +
    ylab("Predicted Observation") +
    ggtitle("Un-shifted Fragments")
}


#' Plot an annealed fragment against the smoothed data.
#'
#' This plot provides a way to assess how good the annealed fragment fits
#' the data.
#'
#' Truncating is useful for excluding the parts of fragments that extend beyond
#' your most recent observation, and therefore do not factor into the
#' loss function.
#'
#' @param data .
#' @param x The datetime column.
#' @param y The column of smoothed (fitted) observations.
#' @param anneal_output Output of `anneal`.
#' @param fragment_obs_col_name The (possibly smoothed) fragments' observations column name.
#' @param k Which fragment number to plot.
#' @param truncate Truncate the plot to `data`'s most recent observation?
#'
#' @export
plot_anneal_fragment <- function(data,
                                 x,
                                 y,
                                 anneal_output,
                                 fragment_obs_col_name,
                                 k = 1,
                                 truncate = TRUE) {
  # Filter to particular k
  k_val = k  # prevent shadowing
  anneal_output$fragments = anneal_output$fragments %>% filter(k == k_val)
  anneal_output$losses = anneal_output$losses %>% filter(k == k_val)

  if (truncate) {
    anneal_output$fragments = anneal_output$fragments %>%
      filter(
        final_datetime < (data %>% pull({{ x }}) %>% max())
      )
  }

  original = anneal_output$fragments %>%
    filter(shift == 0)
  min_losses = anneal_output$losses %>%
    filter(loss == min(loss)) %>%
    pull(shift)
  best = anneal_output$fragments %>%
    filter(shift %in% min_losses)

  ggplot() +
    geom_line(
      data = anneal_output$fragments,
      mapping = aes(
        x = final_datetime,
        y = .data[[fragment_obs_col_name]],
        group = shift,
        color = "Others",
      ),
      alpha = 0.1
    ) +
    geom_line(
      data = original,
      mapping = aes(
        x = final_datetime,
        y = .data[[fragment_obs_col_name]],
        group = shift,
        color = "Original"
      ),
    ) +
    geom_line(
      data = best,
      mapping = aes(
        x = final_datetime,
        y = .data[[fragment_obs_col_name]],
        group = shift,
        color = "Min Loss"
      ),
    ) +
    geom_line(
      data = data,
      mapping = aes(
        x = {{ x }},
        y = {{ y }}
      )
    ) +
    ggtitle(sprintf("Optimal shift: %s", min_losses)) +
    xlab("Datetime") +
    ylab("Predicted Observation") +
    scale_color_manual(
      name = "Fragment",
      values = c(
        "Original" = "red",
        "Min Loss" = "green",
        "Others" = "blue"
      )
    )
}


#' Plot the losses output of anneal.
#'
#' @param anneal_out Output of `anneal`.
#'
#' @export
plot_anneal_losses <- function(anneal_out) {
  # TODO: loss fn attr?
  anneal_out$losses %>%
    mutate(k = as.factor(k)) %>%
    ggplot() +
    geom_point(aes(shift, loss, color = k))
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
