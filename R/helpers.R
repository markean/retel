#' Exponentially tilted empirical likelihood objective function
#'
#' Computes the exponentially tilted empirical likelihood objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param n A single integer.
#' @return A single numeric.
#' @noRd
eval_d_fun <- function(l, g, n) {
  sum(exp(g %*% l)) / n
}

#' Gradient of exponentially tilted empirical likelihood objective function
#'
#' Computes the gradient of the exponentially tilted empirical likelihood
#' objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param n A single integer.
#' @return A single numeric.
#' @noRd
eval_gr_d_fun <- function(l, g, n) {
  colSums(as.vector(exp(g %*% l)) * g) / n
}
