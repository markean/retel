#' Exponentially tilted empirical likelihood
#'
#' Computes exponentially tilted empirical likelihood.
#'
#' @param fn
#'   An estimating function that takes the data `x` and parameter value
#'   `par` as its arguments, returning a numeric matrix. Each row is the return
#'   value from the corresponding row in `x`.
#' @param x
#'   A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation. The number of rows must be
#'   greater than the number of columns.
#' @param par
#'   A numeric vector of parameter values to be tested.
#' @param opts
#'   An optional list with optimization options for [nloptr::nloptr()].
#'   Defaults to `NULL`.
#' @return
#'   A single numeric value representing the log-likelihood ratio. It contains
#'   the optimization results as the attribute `optim`.
#' @references
#'   Schennach, SM (2005).
#'   "Bayesian Exponentially Tilted Empirical Likelihood."
#'   \emph{Biometrika}, 92, 31--46.
#' @examples
#' # Generate data
#' set.seed(63456)
#' x <- rnorm(100)
#'
#' # Define an estimating function (ex. mean)
#' fn <- function(x, par) {
#'   x - par
#' }
#'
#' # Set parameter value
#' par <- 0
#'
#' # Call the etel function
#' etel(fn, x, par)
#' @export
etel <- function(fn, x, par, opts = NULL) {
  assert_function(fn, args = c("x", "par"), ordered = TRUE, nargs = 2L)

  if (isFALSE(is.matrix(x))) {
    x <- as.matrix(x, rownames.force = TRUE)
  }
  assert_matrix(x,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 2L,
    min.cols = 1L
  )
  n <- nrow(x)

  assert_numeric(par,
    finite = TRUE, any.missing = FALSE, all.missing = FALSE, min.len = 1L,
    typed.missing = TRUE
  )

  g <- fn(x, par)
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g, rownames.force = TRUE)
  }
  assert_matrix(g,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 2L,
    min.cols = 1L, nrows = n
  )
  p <- ncol(g)
  stopifnot(
    "`fn(x, par)` must produce a numeric matrix with column rank." =
      (isTRUE(n > p && rankMatrix(g) == p))
  )

  if (isTRUE(is.null(opts))) {
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
  } else {
    assert_list(opts)
  }

  optim <- nloptr(
    x0 = rep(0, p), eval_f = eval_d_fn, eval_grad_f = eval_gr_d_fn, opts = opts,
    g = g, n = n, tau = 0
  )
  lambda <- optim$solution

  out <- as.numeric(lambda %*% colSums(g)) - n * log(eval_d_fn(lambda, g, n, 0))
  attributes(out) <- list(optim = optim)
  out
}
