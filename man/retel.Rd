% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retel.R
\name{retel}
\alias{retel}
\title{Regularized exponentially tilted empirical likelihood}
\usage{
retel(fn, x, par, mu, Sigma, tau, type = "full", opts = NULL)
}
\arguments{
\item{fn}{An estimating function that takes the data \code{x} and parameter value
\code{par} as its arguments, returning a numeric matrix. Each row is the return
value from the corresponding row in \code{x}.}

\item{x}{A numeric matrix, or an object that can be coerced to a numeric
matrix. Each row corresponds to an observation. The number of rows must be
greater than the number of columns.}

\item{par}{A numeric vector of parameter values to be tested.}

\item{mu}{A numeric vector of parameters for regularization. See 'Details' for more
information.}

\item{Sigma}{A numeric matrix, or an object that can be coerced to a numeric matrix,
of parameters for regularization. See 'Details' for more information.}

\item{tau}{A single numeric parameter for regularization. See 'Details' for more
information.}

\item{type}{A single character indicating the type of regularized exponentially tilted
empirical likelihood. It must be either \code{"full"} or \code{"reduced"}. Defaults
to \code{"full"}.  See 'Details' for more information.}

\item{opts}{An optional list with optimization options for \code{\link[nloptr:nloptr]{nloptr::nloptr()}}.
Defaults to \code{NULL}.}
}
\value{
A single numeric value representing the log-likelihood ratio. It contains
the optimization results as the attribute \code{optim}.
}
\description{
Computes regularized exponentially tilted empirical likelihood.
}
\details{
Let \eqn{\{\bm{X}_i\}_{i = 1}^n} denote independent \eqn{d_x}-dimensional
observations from a complete probability space
\eqn{{(\mathcal{X}, \mathcal{F}, P)}} satisfying the moment condition:
\deqn{\textnormal{E}_P[\bm{g}(\bm{X}_i, \bm{\theta})] = \bm{0},}
where \eqn{{\bm{g}}:
  {\mathbb{R}^{d_x} \times \Theta} \mapsto {\mathbb{R}^p}} is an estimating
function with the true parameter value
\eqn{{\bm{\theta}_0} \in {\Theta} \subset \mathbb{R}^p}.

For a given parameter value \eqn{\bm{\theta}}, regularized exponentially
tilted empirical likelihood solves the following optimization problem:
\deqn{
    \min_{\bm{\lambda} \in \mathbb{R}^p}
      \left\{
        d_n\left(\bm{\theta}, \bm{\lambda}\right) +
        p_n\left(\bm{\theta}, \bm{\lambda}\right)
      \right\},
  }
where
\deqn{
    d_n\left(\bm{\theta}, \bm{\lambda}\right) =
    \frac{1}{n + \tau_n}
    \sum_{i = 1}^n \exp
      \left(
        \bm{\lambda}^\top \bm{g}\left(\bm{X}_i, \bm{\theta}\right)
      \right)
  }
and
\deqn{
    p_n\left(\bm{\theta}, \bm{\lambda}\right) =
    \frac{\tau_n}{n + \tau_n}
    \exp
      \left(
        \bm{\lambda}^\top\bm{\mu}_{n, \bm{\theta}} +
        \frac{1}{2}
        \bm{\lambda}^\top\bm{\Sigma}_{n, \bm{\theta}}\bm{\lambda}
      \right).
  }
Here, \eqn{{\tau_n} > {0}}, \eqn{\bm{\mu}_{n, \bm{\theta}}},
\eqn{\bm{\Sigma}_{n, \bm{\theta}}} are all tuning parameters that control
the strength of \eqn{{p_n(\bm{\theta}, \bm{\lambda})}} as a penalty.

Once we have determined the solution \eqn{{\bm{\lambda}_{RET}}}, we define
the likelihood ratio function as follows:
\deqn{
    R_{RET}\left(\bm{\theta}\right) =
    \left(
      \frac{n + \tau_n}{\tau_n}p_c\left(\bm{\theta}\right)\right)
      \prod_{i = 1}^n \left(n + \tau_n\right)p_i\left(\bm{\theta}
    \right),
  }
where
\deqn{
    p_i\left(\bm{\theta}\right) =
    \frac{\exp
      \left(
        {\bm{\lambda}_{RET}}^\top\bm{g}\left(\bm{X}_i, \bm{\theta}\right)
      \right)
    }{c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
    } \quad \left(i = 1, \dots, n\right),\quad
    p_c\left(\bm{\theta}\right) =
    \frac{p_n
      \left(\bm{\theta}, \bm{\lambda}_{RET}\right)
    }{c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
    },
  }
and
\eqn{
    c_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right) =
    d_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right) +
    p_n\left(\bm{\theta}, \bm{\lambda}_{RET}\right)
  }. The reduced version of the likelihood ratio function is defined as:
\deqn{
    \widetilde{R}_{RET}\left(\bm{\theta}\right) =
    \prod_{i = 1}^n \left(n + \tau_n\right)p_i\left(\bm{\theta}\right).
  }

See the references below for more details on derivation, interpretation,
and properties.
}
\examples{
# Generate data
set.seed(63456)
x <- rnorm(100)

# Define an estimating function (ex. mean)
fn <- function(x, par) {
  x - par
}

# Set parameter value
par <- 0

# Set regularization parameters
mu <- 0
Sigma <- 1
tau <- 1

# Call the retel function
retel(fn, x, par, mu, Sigma, tau)
}
\references{
Kim E, MacEachern SN, Peruggia M (2023).
"Regularized Exponentially Tilted Empirical Likelihood for Bayesian
Inference." \doi{10.48550/arXiv.2312.17015}.
}
