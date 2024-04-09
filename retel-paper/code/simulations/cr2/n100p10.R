## 1. Load packages
library(ggplot2)
library(melt)
library(mvtnorm)
library(nloptr)
library(retel)


## 2. Parameters
n <- 100
p <- 10
an <- log(n) / 2
par <- rep(0, p)
cv <- qchisq(0.95, df = p)


## 3. Constants
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
S <- 1e+03L


## 4. Functions
f <- function(x, par) {
  x - par
}
obj_aetel <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  mean(exp(g %*% l))
}
gr_obj_aetel <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  colMeans(as.vector(exp(g %*% l)) * g)
}
aetel <- function(x, theta) {
  g <- c(x - theta, -log(n) / 2 * mean(x - theta))
  out <- nloptr(
    x0 = 0, eval_f = obj_aetel, eval_grad_f = gr_obj_aetel, opts = opts, g = g
  )
  lambda <- out$solution
  -(n + 1) * log(obj_aetel(lambda, g)) + lambda * sum(g)
}


## 5. Simulation
set.seed(6325203)
cat("< Simulation replications =", S, ">\n")
cat("< Sample size =", n, ">\n")
cat("< Parameter dimension =", p, ">\n")
res_retel_f <- vector("numeric", length = S)
res_retel_r <- vector("numeric", length = S)
res_aetel <- vector("numeric", length = S)
res_ael <- vector("numeric", length = S)
res_bael <- vector("numeric", length = S)
cov_retel_f <- vector("logical", length = S)
cov_retel_r <- vector("logical", length = S)
cov_aetel <- vector("logical", length = S)
cov_ael <- vector("logical", length = S)
cov_bael <- vector("logical", length = S)
for (i in seq_len(S)) {
  # Data
  x <- rmvnorm(n, sigma = diag(p))

  # RETEL_f
  res_retel_f[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "full"
  )
  cov_retel_f[i] <- res_retel_f[i] < cv

  # RETEL_r
  res_retel_r[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "reduced"
  )
  cov_retel_r[i] <- res_retel_r[i] < cv

  # AETEL
  res_aetel[i] <- -2 * aetel(x, par)
  cov_aetel[i] <- res_aetel[i] < cv

  # AEL
  p_ael <- par - an * (colMeans(x) - par)
  res_ael[i] <- el_mean(rbind(x, p_ael), par = par) |>
    chisq()
  cov_ael[i] <- res_ael[i] < cv

  # BAEL
  v <- colMeans(x) - par
  u <- v / norm(v, type = "2")
  p1_bael <- par - 1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  p2_bael <- 2 * colMeans(x) - par +
    1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  res_bael[i] <- el_mean(rbind(x, p1_bael, p2_bael), par = par) |>
    chisq()
  cov_bael[i] <- res_bael[i] < cv
}


## 5. Results
mean(cov_retel_f)
mean(cov_retel_r)
mean(cov_aetel)
mean(cov_ael)
mean(cov_bael)
