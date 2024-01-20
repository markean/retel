## 1. Load packages
options(warn = -1)
options(scipen = 999)
suppressMessages(library(actuar))
library(coda)
suppressMessages(library(here))
suppressMessages(i_am("code/application/etel.R"))
library(mvtnorm)
library(retel)


## 2. Load data
data("income", package = "retel")
raw_income <- read.csv(paste0(dirname(getwd()), "/data-raw/income.csv"))


## 3. Constants
# Model
fit <- lm(mi_1989 ~ -1 + mi_1979 + ami, income)
b0 <- coef(fit)
x <- model.matrix(fit)
y <- as.numeric(model.response(fit$model))
m <- nrow(x)
p <- ncol(x)
# MCMC parameters
B <- 250000L


## 4. Functions
f <- function(x, par) {
  cbind(x - par, (x - par)^2L - 1)
}
# Log prior density functions
log_pd_theta <- function(theta, b, s2) {
  dmvnorm(theta, mean = x %*% b, sigma = s2 * diag(m), log = TRUE)
}
log_pd_beta <- function(b, s2) {
  dmvnorm(b, mean = b0, sigma = 0.1 * s2 * solve(crossprod(x)), log = TRUE)
}
log_pd_s2 <- function(s2) {
  dinvgamma(s2, shape = 4, scale = 1, log = TRUE)
}
log_pd <- function(theta, b, s2) {
  log_pd_theta(theta, b, s2) + log_pd_beta(b, s2) - log(s2)
}
log_posterior_unnormalized <- function(theta, b, s2) {
  log_pd(theta, b, s2) + etel(f, y, theta)
}
# MCMC
mcmc_fn <- function(B) {
  theta_sample <- matrix(nrow = B, ncol = m)
  theta_sample[1L, ] <- rnorm(m, mean = mean(y), sd = 0.1)
  beta_sample <- matrix(nrow = B, ncol = p)
  beta_sample[1L, ] <- rnorm(p, mean = b0, sd = 0.1)
  s2_sample <- vector("numeric", length = B)
  s2_sample[1L] <- rinvgamma(1L, shape = 4, scale = 1)
  acceptace <- vector("logical", length = B)
  acceptace[1L] <- FALSE
  for (i in seq_len(B)[-1]) {
    # Sample proposal value
    s2_p <- rnorm(1L, mean = s2_sample[i - 1L], sd = 0.06)
    if (s2_p < 0) {
      s2_p <- 0
    }
    b_p <- rmvnorm(1L,
      mean = beta_sample[i - 1L, ], sigma = 0.01 * diag(p)
    ) |>
      as.vector()
    theta_p <- rmvnorm(1L,
      mean = theta_sample[i - 1L, ], sigma = 0.06 * diag(m)
    ) |>
      as.vector()
    # Compute log ratio of unnormailzed posterior densities
    logr <- log_posterior_unnormalized(theta_p, b_p, s2_p) -
      log_posterior_unnormalized(
        theta_sample[i - 1L, ], beta_sample[i - 1L, ],
        s2_sample[i - 1L]
      )
    # Sample uniform random variable
    u <- runif(1L)
    # Accept or reject
    if (isTRUE(log(u) < logr)) {
      theta_sample[i, ] <- theta_p
      beta_sample[i, ] <- b_p
      s2_sample[i] <- s2_p
      acceptace[i] <- TRUE
    } else {
      theta_sample[i, ] <- theta_sample[i - 1L, ]
      beta_sample[i, ] <- beta_sample[i - 1L, ]
      s2_sample[i] <- s2_sample[i - 1L]
      acceptace[i] <- FALSE
    }
  }
  list(
    theta = theta_sample, beta = beta_sample, s2 = s2_sample, rate = acceptace
  )
}
# Metrics
aad <- function(y, e) {
  mean(abs(y - e))
}
aard <- function(y, e) {
  mean(abs((y - e) / y))
}
asd <- function(y, e) {
  mean((y - e)^2L)
}
asrd <- function(y, e) {
  mean(((y - e) / y)^2L)
}


## 5. Simulations
set.seed(2536345)
cat("< Method = ETEL>\n")
cat("< MCMC runs =", 4L * B, ">\n")
c1 <- mcmc_fn(B)
c2 <- mcmc_fn(B)
c3 <- mcmc_fn(B)
c4 <- mcmc_fn(B)
theta_c1 <- mcmc(c1$theta)
theta_c2 <- mcmc(c2$theta)
theta_c3 <- mcmc(c3$theta)
theta_c4 <- mcmc(c4$theta)
beta_c1 <- mcmc(c1$beta)
beta_c2 <- mcmc(c2$beta)
beta_c3 <- mcmc(c3$beta)
beta_c4 <- mcmc(c4$beta)
s2_c1 <- mcmc(c1$s2)
s2_c2 <- mcmc(c2$s2)
s2_c3 <- mcmc(c3$s2)
s2_c4 <- mcmc(c4$s2)


## 6. Results
theta <- rbind(c1$theta, c2$theta, c3$theta, c4$theta)
beta <- rbind(c1$beta, c2$beta, c3$beta, c4$beta)
s2 <- c(c1$s2, c2$s2, c3$s2, c4$s2)
accept <- c(c1$rate, c2$rate, c4$rate, c4$rate)
# Potential scale reduction factors
cat("\ntheta psrf:")
cat("\n")
gelman.diag(mcmc.list(theta_c1, theta_c2, theta_c3, theta_c4))$psrf[, 1L]
cat(
  "theta mpsrf:",
  gelman.diag(mcmc.list(theta_c1, theta_c2, theta_c3, theta_c4))$mpsrf
)
cat(
  "\nbeta psrf:",
  gelman.diag(mcmc.list(beta_c1, beta_c2, beta_c3, beta_c4))$psrf[, 1L]
)
cat(
  "\nbeta mpsrf:",
  gelman.diag(mcmc.list(beta_c1, beta_c2, beta_c3, beta_c4))$mpsrf
)
cat(
  "\ns2 psrf:",
  gelman.diag(mcmc.list(s2_c1, s2_c2, s2_c3, s2_c4))$psrf[1L]
)
# Acceptance rate
cat("\nAcceptance rate: ", mean(mean(c(c1$rate, c2$rate, c4$rate, c4$rate))))
cat("\n\n")
# Summary of beta and s2
cat("< beta & s2 summary >\n")
cat("beta1 mean:   ", mean(beta[, 1L]))
cat("\nbeta1 95% ci: ", quantile(beta[, 1L], c(0.025, 0.975)))
cat("\nbeta2 mean:   ", mean(beta[, 2L]))
cat("\nbeta2 95% ci: ", quantile(beta[, 2L], c(0.025, 0.975)))
cat("\ns2 mean:      ", mean(s2))
cat("\ns2 95% ci:    ", quantile(s2, c(0.025, 0.975)))
# Summary of theta
theta_median <- apply(theta, MARGIN = 2L, FUN = median)
theta_ci <- t(apply(theta, MARGIN = 2L, FUN = quantile, c(0.025, 0.975)))
theta_al <- mean(theta_ci[, 2L] - theta_ci[, 1L])
cat("\n\n< theta summary >\n")
cat("theta 95% ci al: ", theta_al)
cat("\ntheta aad:       ", aad(y, theta_median))
cat("\ntheta aard:      ", aard(y, theta_median))
cat("\ntheta asd:       ", asd(y, theta_median))
cat("\ntheta asrd:      ", asrd(y, theta_median))
# Original scale
y_raw <- raw_income$mi_1989
theta_os <- t(apply(theta,
  MARGIN = 1L, FUN = function(x) x * sd(y_raw) + mean(y_raw)
))
theta_os_median <- apply(theta_os, MARGIN = 2L, FUN = median)
theta_os_ci <- t(apply(theta_os, MARGIN = 2L, FUN = quantile, c(0.025, 0.975)))
theta_os_al <- mean(theta_os_ci[, 2L] - theta_os_ci[, 1L])
cat("\n\n< theta_os summary >\n")
cat("theta_os 95% ci al: ", theta_os_al)
cat("\ntheta_os aad:       ", aad(y_raw, theta_os_median))
cat("\ntheta_os aard:      ", aard(y_raw, theta_os_median))
cat("\ntheta_os asd:       ", asd(y_raw, theta_os_median))
cat("\ntheta_os asrd:      ", asrd(y_raw, theta_os_median))
cat("\n")
