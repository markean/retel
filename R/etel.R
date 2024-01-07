#' Exponentially tilted empirical likelihood
#'
#' Computes exponentially tilted empirical likelihood.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating function.
#'   The number of rows must be greater than the number of columns.
#' @return An object of class?.
#' @references Schennach, SM (2005).
#'   "Bayesian exponentially tilted empirical likelihood"
#'   \emph{Biometrika}, 92, 31--46.
#' @export
etel <- function(g) {
  mm <- as.matrix(g, rownames.force = TRUE)
  n <- nrow(mm)
  p <- ncol(mm)
  stopifnot(
    "`g` must have at least two observations." = (n >= 2L),
    "`g` must be a finite numeric matrix." =
      (isTRUE(is.numeric(mm) && all(is.finite(mm)))),
    "`g` must have full column rank." = (isTRUE(n > p && rankMatrix(mm) == p))
  )
  obj <- function(l, g) {
    sum(exp(g %*% l)) / n
  }
  gr_obj <- function(l, g) {
    colMeans(as.vector(exp(g %*% l)) * g)
  }
  opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
  out <- nloptr(
    x0 = rep(0, NCOL(g)), eval_f = obj, eval_grad_f = gr_obj, opts = opts,
    g = mm
  )
  print(out$status)
  print(out$message)
  lambda <- out$solution
  -n * log(obj(lambda, mm)) + as.numeric(lambda %*% colSums(mm))
}
