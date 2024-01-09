#' Exponentially tilted empirical likelihood
#'
#' Computes exponentially tilted empirical likelihood.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating function.
#'   The number of rows must be greater than the number of columns.
#' @param opts A list with optimization options for [nloptr()].
#' @return A single numeric of the log-likelihood ratio.
#' @references Schennach, SM (2005).
#'   "Bayesian exponentially tilted empirical likelihood"
#'   \emph{Biometrika}, 92, 31--46.
#' @examples
#' set.seed(63456)
#' x <- rnorm(100)
#' par <- 0
#' g <- x - par
#' etel(g)
#' @export
etel <- function(g, opts) {
  mm <- as.matrix(g, rownames.force = TRUE)
  n <- nrow(mm)
  p <- ncol(mm)
  stopifnot(
    "`g` must have at least two observations." = (n >= 2L),
    "`g` must be a finite numeric matrix." =
      (isTRUE(is.numeric(mm) && all(is.finite(mm)))),
    "`g` must have full column rank." = (isTRUE(n > p && rankMatrix(mm) == p))
  )
  if (missing(opts)) {
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)
  }
  optim <- nloptr(
    x0 = rep(0, NCOL(g)), eval_f = eval_d_fun, eval_grad_f = eval_gr_d_fun,
    opts = opts, g = mm, n = n
  )
  lambda <- optim$solution
  out <-
    -n * log(eval_d_fun(lambda, mm, n)) + as.numeric(lambda %*% colSums(mm))
  attributes(out) <- list(optim = optim)
  out
}
