#' Regularized exponentially tilted empirical likelihood
#'
#' Computes regularized exponentially tilted empirical likelihood.
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
#' @param mu
#'   A numeric vector of parameters for regularization. See 'Details' for more
#'   information.
#' @param Sigma
#'   A numeric matrix, or an object that can be coerced to a numeric matrix,
#'   of parameters for regularization. See 'Details' for more information.
#' @param tau
#'   A single numeric parameter for regularization. See 'Details' for more
#'   information.
#' @param type
#'   A single character indicating the type of regularized exponentially tilted
#'   empirical likelihood. It must be either `"full"` or `"reduced"`. Defaults
#'   to `"full"`.  See 'Details' for more information.
#' @param opts
#'   An optional list with optimization options for [nloptr()].
#'   Defaults to `NULL`.
#' @details
#'   Let \eqn{\{\bm{X}_i\}_{i = 1}^n} denote independent \eqn{d_x}-dimensional
#'   observations from a complete probability space
#'   \eqn{{(\mathcal{X}, \mathcal{F}, P)}} satisfying the moment condition:
#'   \deqn{\textnormal{E}_P[\bm{g}(\bm{X}_i, \bm{\theta})] = \bm{0},}
#'   where \eqn{{\bm{g}}:
#'   {\mathbb{R}^{d_x} \times \Theta} \mapsto {\mathbb{R}^p}} is an estimating
#'   function with the true parameter value
#'   \eqn{{\bm{\theta}_0} \in {\Theta} \subset \mathbb{R}^p}.
#'
#'   For a given parameter value \eqn{\bm{\theta}}, regularized exponentially
#'   tilted empirical likelihood solves the following optimization problem:
#'   \deqn{
#'     \min_{\bm{\lambda} \in \mathbb{R}^p}
#'       \left\{
#'         d_n\left(\bm{\theta}, \bm{\lambda}\right) +
#'         p_n\left(\bm{\theta}, \bm{\lambda}\right)
#'       \right\},
#'   }
#'   where
#'   \deqn{
#'     d_n\left(\bm{\theta}, \bm{\lambda}\right) =
#'     \frac{1}{n + \tau_n}
#'     \sum_{i = 1}^n \exp
#'       \left(
#'         \bm{\lambda}^\top \bm{g}\left(\bm{X}_i, \bm{\theta}\right)
#'       \right)
#'   }
#'   and
#'   \deqn{
#'     p_n\left(\bm{\theta}, \bm{\lambda}\right) =
#'     \frac{\tau_n}{n + \tau_n}
#'     \exp
#'       \left(
#'         \bm{\lambda}^\top\bm{\mu}_{n, \bm{\theta}} +
#'         \frac{1}{2}
#'         \bm{\lambda}^\top\bm{\Sigma}_{n, \bm{\theta}}\bm{\lambda}
#'       \right).
#'   }
#'   Here, \eqn{{\tau_n} > {0}}, \eqn{\bm{\mu}_{n, \bm{\theta}}},
#'   \eqn{\bm{\Sigma}_{n, \bm{\theta}}} are all tuning parameters that control
#'   the strength of \eqn{{p_n(\bm{\theta}, \bm{\lambda})}} as a penalty.
#'
#'   Once we have determined the solution \eqn{{\bm{\lambda}_{RET}}}, we define
#'   the likelihood ratio function as follows:
#'   \deqn{
#'     R_{RET}\left(\bm{\theta}\right) =
#'     \left(
#'       \frac{n + \tau_n}{\tau_n}p_c\left(\bm{\theta}\right)\right)
#'       \prod_{i = 1}^n \left(n + \tau_n\right)p_i\left(\bm{\theta}
#'     \right),
#'   }
#'   where
#'   \deqn{
#'     p_i\left(\bm{\theta}\right) =
#'     \frac{\exp
#'       \left(
#'         {\bm{\lambda}_{RET}}^\top\bm{g}\left(\bm{X}_i, \bm{\theta}\right)
#'       \right)
#'     }{c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
#'     } \quad \left(i = 1, \dots, n\right),\quad
#'     p_c\left(\bm{\theta}\right) =
#'     \frac{p_n
#'       \left(\bm{\theta}, \bm{\lambda}_{RET}\right)
#'     }{c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
#'     },
#'   }
#'   and
#'   \eqn{
#'     c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right) =
#'     d_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right) +
#'     p_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
#'   }. The reduced version of the likelihood ratio function is defined as:
#'   \deqn{
#'     \widetilde{R}_{RET}\left(\bm{\theta}\right) =
#'     \prod_{i = 1}^n \left(n + \tau_n\right)p_i\left(\bm{\theta}\right).
#'   }
#'
#'   See the references below for more details on derivation, interpretation,
#'   and properties.
#' @return
#'   A single numeric value representing the log-likelihood ratio. It contains
#'   the optimization results as the attribute `optim`.
#' @references
#'   Kim E, MacEachern SN, Peruggia M (2023).
#'   "Regularized Exponentially Tilted Empirical Likelihood for Bayesian
#'   Inference." \doi{10.48550/arXiv.2312.17015}.
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
#' # Set regularization parameters
#' mu <- 0
#' Sigma <- 1
#' tau <- 1
#'
#' # Call the retel function
#' retel(fn, x, par, mu, Sigma, tau)
#' @export
retel <- function(fn, x, par, mu, Sigma, tau, type = "full", opts = NULL) {
  assert_function(fn, args = c("x", "par"), ordered = TRUE, nargs = 2L)

  if (isFALSE(is.matrix(x))) {
    x <- as.matrix(x, rownames.force = TRUE)
  }
  assert_matrix(x,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 1L,
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
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 1L,
    min.cols = 1L, nrows = n
  )
  p <- ncol(g)
  stopifnot(
    "`fn(x, par)` must produce a numeric matrix with column rank." =
      (isTRUE(n >= p && rankMatrix(g) == p))
  )

  assert_numeric(mu,
    finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = p,
    typed.missing = TRUE
  )

  if (isFALSE(is.matrix(Sigma))) {
    Sigma <- as.matrix(Sigma, rownames.force = TRUE)
  }
  assert_matrix(Sigma,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, nrows = p,
    ncols = p
  )
  stopifnot(
    "`Sigma` must be a positive definite matrix." =
      (isTRUE(is.positive.definite(Sigma)))
  )

  assert_number(tau, lower = 1, finite = TRUE)

  assert_choice(type, c("full", "reduced"))

  if (isTRUE(is.null(opts))) {
    opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
  } else {
    assert_list(opts)
  }

  optim <- nloptr(
    x0 = rep(0, p), eval_f = eval_obj_fn, eval_grad_f = eval_gr_obj_fn,
    opts = opts, g = g, mu = mu, Sigma = Sigma, n = n, tau = tau
  )
  lambda <- optim$solution

  if (isTRUE(type == "full")) {
    out <- as.numeric(lambda %*% colSums(g)) + as.numeric(lambda %*% mu) +
      colSums(lambda * (Sigma %*% lambda)) / 2 -
      (n + 1) * log(eval_obj_fn(lambda, g, mu, Sigma, n, tau))
  } else {
    out <- as.numeric(lambda %*% colSums(g)) -
      n * log(eval_obj_fn(lambda, g, mu, Sigma, n, tau))
  }
  attributes(out) <- list(optim = optim)
  out
}
