#' Regularized exponentially tilted empirical likelihood
#'
#' Computes regularized exponentially tilted empirical likelihood.
#'
#' @param fn An estimating function that takes the data `x` and parameter value
#'   `par` as its arguments, returning a numeric matrix. Each row is the return
#'   value from the corresponding row in `x`.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation. The number of rows must be
#'   greater than the number of columns.
#' @param par A numeric vector of parameter values to be tested. The length of
#'   the vector must be the same as the number of columns in `x`.
#' @param mu A numeric matrix.
#' @param Sigma A numeric matrix.
#' @param tau A single numeric.
#' @param opts A list with optimization options for [nloptr()].
#' @return A single numeric of the log-likelihood ratio.
#' @references Kim E, MacEachern SN, Peruggia M (2023).
#'   "Regularized Exponentially Tilted Empirical Likelihood for Bayesian
#'   Inference."
#'   \doi{10.48550/arXiv.2312.17015}.
#' @examples
#' set.seed(63456)
#' f <- function(x, par) {
#'   x - par
#' }
#' x <- rnorm(100)
#' par <- 0
#' mu <- 0
#' Sigma <- matrix(rnorm(1), nrow = 1)
#' tau <- 1
#' retel(f, x, par, mu, Sigma, tau)
#' @export
retel <- function(fn, x, par, mu, Sigma, tau, opts) {
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
    x0 = rep(0, ncol(g)), eval_f = eval_obj_fn, eval_grad_f = eval_gr_obj_fn,
    opts = opts, g = g, mu = mu, Sigma = Sigma, tau = tau, n = n
  )
  lambda <- optim$solution
  out <- as.numeric(lambda %*% colSums(g)) + as.numeric(lambda %*% mu) +
    colSums(mu * (Sigma %*% mu)) / 2 -
    (n + 1) * log(n) + (n + 1) * log(n + tau) -
    (n + 1) * log(eval_obj_fn(lambda, g, mu, Sigma, tau, n))
  attributes(out) <- list(optim = optim)
  out
}
