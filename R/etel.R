#' Exponentially tilted empirical likelihood
#'
#' Computes exponentially tilted empirical likelihood.
#'
#' @param fn An estimating function that takes the data `x` and parameter value
#'   `par` as its arguments, returning a numeric matrix. Each row is the return
#'   value from the corresponding row in `x`.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation. The number of rows must be
#'   greater than the number of columns.
#' @param par A numeric vector of parameter values to be tested. The length of
#'   the vector must be the same as the number of columns in `x`.
#' @param opts A list with optimization options for [nloptr()].
#' @return A single numeric of the log-likelihood ratio.
#' @references Schennach, SM (2005).
#'   "Bayesian exponentially tilted empirical likelihood"
#'   \emph{Biometrika}, 92, 31--46.
#' @examples
#' set.seed(63456)
#' f <- function(x, par) {
#'   x - par
#' }
#' x <- rnorm(100)
#' par <- 0
#' etel(f, x, par)
#' @export
etel <- function(fn, x, par, opts) {
  x <- validate_x(x)
  par <- validate_par(par)
  g <- fn(x, par)
  g <- as.matrix(g, rownames.force = TRUE)
  n <- nrow(g)
  p <- ncol(g)
  stopifnot(
    "`g` must have at least two observations." = (n >= 2L),
    "`g` must be a finite numeric matrix." =
      (isTRUE(is.numeric(g) && all(is.finite(g)))),
    "`g` must have full column rank." = (isTRUE(n > p && rankMatrix(g) == p))
  )
  if (missing(opts)) {
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)
  }
  optim <- nloptr(
    x0 = rep(0, ncol(g)), eval_f = eval_d_fn, eval_grad_f = eval_gr_d_fn,
    opts = opts, g = g, n = n, tau = 0
  )
  lambda <- optim$solution
  out <- as.numeric(lambda %*% colSums(g)) - n * log(eval_d_fn(lambda, g, n, 0))
  attributes(out) <- list(optim = optim)
  out
}
