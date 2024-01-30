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


#' Calculate the root mean squared error of two sequences, weighted by slope.
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
#' extremeness of summer/winter.
#'
#' This is not always the best option, however, and results should be
#' inspected.
#'
#' @param a Sequence 1
#' @param b Sequence 2
#' @returns Slope-weighted RMSE
#'
#' @export
slope_weighted_rmse <- function(a, b) {
  b_diff = difference(b)

  sqrt(mean(((a - b) * b_diff) ^ 2, na.rm = TRUE))
}


#' Calculate the root mean squared error of two sequences, weighted by recency.
#'
#' Weights are based on the length(a) linear sequence over [0,1].
#' At time t, the formulation is: ((a_t - b_t) ^ 2) * wt_t.
#'
#' Intended to be passed as the `loss_fn` arg of `anneal`.
#'
#' In the context of this package, this punishes strong deviations near the
#' head of your time series. For instance, a sudden flattening in the last
#' few observations may be important.
#'
#' This is not always the best option, however, and results should be
#' inspected.
#'
#' @param a Sequence 1
#' @param b Sequence 2
#' @returns Recency-weighted RMSE
#'
#' @export
recency_weighted_rmse <- function(a, b) {
  distance_wts = seq(0, 1, length.out = length(a))
  sqrt(mean(((a - b) ^ 2) * distance_wts, na.rm = TRUE))
}
