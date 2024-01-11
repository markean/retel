#' Exponentially tilted empirical likelihood objective function
#'
#' Computes the exponentially tilted empirical likelihood objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A single numeric.
#' @noRd
eval_d_fn <- function(l, g, n, tau) {
  sum(exp(g %*% l)) / (n + tau)
}

#' Exponentially tilted empirical likelihood objective function
#'
#' Computes the exponentially tilted empirical likelihood objective function.
#'
#' @param l A numeric vector.
#' @param mu A numeric matrix.
#' @param Sigma A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A single numeric.
#' @noRd
eval_p_fn <- function(l, mu, Sigma, n, tau) {
  tau * exp(as.numeric(l %*% mu) + colSums(l * (Sigma %*% l)) / 2) / (n + tau)
}

#' Exponentially tilted empirical likelihood objective function
#'
#' Computes the exponentially tilted empirical likelihood objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param mu A numeric matrix.
#' @param Sigma A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A single numeric.
#' @noRd
eval_obj_fn <- function(l, g, mu, Sigma, n, tau) {
  eval_d_fn(l, g, n, tau) + eval_p_fn(l, mu, Sigma, n, tau)
}

#' Gradient of exponentially tilted empirical likelihood objective function
#'
#' Computes the gradient of the exponentially tilted empirical likelihood
#' objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A numeric vector.
#' @noRd
eval_gr_d_fn <- function(l, g, n, tau) {
  colSums(as.numeric(exp(g %*% l)) * g) / (n + tau)
}

#' Gradient of exponentially tilted empirical likelihood objective function
#'
#' Computes the gradient of the exponentially tilted empirical likelihood
#' objective function.
#'
#' @param l A numeric vector.
#' @param mu A numeric matrix.
#' @param Sigma A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A numeric vector.
#' @noRd
eval_gr_p_fn <- function(l, mu, Sigma, n, tau) {
  eval_p_fn(l, mu, Sigma, n, tau) * (mu + as.numeric(Sigma %*% l))
}

#' Gradient of exponentially tilted empirical likelihood objective function
#'
#' Computes the gradient of the exponentially tilted empirical likelihood
#' objective function.
#'
#' @param l A numeric vector.
#' @param g A numeric matrix.
#' @param mu A numeric matrix.
#' @param Sigma A numeric matrix.
#' @param n A single integer.
#' @param tau A single numeric.
#' @return A numeric vector.
#' @noRd
eval_gr_obj_fn <- function(l, g, mu, Sigma, n, tau) {
  eval_gr_d_fn(l, g, n, tau) + eval_gr_p_fn(l, mu, Sigma, n, tau)
}

#' Validate `x`
#'
#' Validate `x` in [etel()] and its extensions.
#'
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix.
#' @return A numeric matrix.
#' @noRd
validate_x <- function(x) {
  x <- as.matrix(x, rownames.force = TRUE)
  stopifnot(
    "`x` must have at least two observations." = (nrow(x) >= 2L),
    "`x` must must have larger number of rows than columns." =
      nrow(x) > ncol(x),
    "`x` must be a finite numeric matrix." =
      isTRUE(is.numeric(x) && all(is.finite(x))),
    "`x` must have full column rank." = rankMatrix(x) == ncol(x)
  )
  x
}

#' Validate `par`
#'
#' Validate `par` in [etel()] and its extensions.
#'
#' @param par A numeric vector.
#' @return A numeric vector.
#' @noRd
validate_par <- function(par) {
  stopifnot(
    "`par` must be a finite numeric vector." =
      (isTRUE(is.numeric(par) && all(is.finite(par))))
  )
  par
}
