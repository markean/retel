## 1. Load packages
options(warn = -1)
options(scipen = 999)
library(coda)
suppressMessages(library(doParallel)) # parallel backend
suppressMessages(library(doRNG)) # reproducible parallel loop
library(melt)
library(mvtnorm)


## 2. Parameters
# Sample size for each group
n <- 10L


## 3. Constants
# Standard deviation of prior
sd_mu <- 10
# Simulation replications
S <- 1e+03L
# MCMC parameters
B <- 5000L
# Optimization
ctrl <- el_control(th = 1000, maxit_l = 50L, nthreads = 1L)


## 4. Functions
# Prior and posterior computation
log_prior <- function(mu, theta) {
  dnorm(mu, mean = 0, sd = sd_mu, log = TRUE) -
    sum(log(pi * (1 + (theta - mu)^2L)))
}
loglr <- function(theta, x1, x2) {
  fit1 <- el_mean(x1, par = theta[1L], control = ctrl)
  fit2 <- el_mean(x2, par = theta[2L], control = ctrl)
  loglr1 <- if (conv(fit1)) logLR(fit1) else -10000
  loglr2 <- if (conv(fit2)) logLR(fit2) else -10000
  loglr1 + loglr2
}
log_posterior_unnormalized <- function(mu_p, theta, x1, x2) {
  log_prior(mu_p, theta) + loglr(theta, x1, x2)
}


## 5. Simulations
cat("< Method = EL>\n")
cat("< Simulation replications =", S, ">\n")
cat("< MCMC runs =", 2L * B, ">\n")
cat("< n =", n, ">\n")
set.seed(1254272)
cl <- makeCluster(24L)
registerDoParallel(cl)
result <- foreach(
  i = icount(S), .combine = "rbind", .inorder = FALSE,
  .packages = c("coda", "melt", "mvtnorm")
) %dorng% {
  # Generate data
  x1 <- rnorm(n, mean = -3, sd = 1)
  x2 <- rnorm(n, mean = 3, sd = 1)

  # MCMC chain 1
  mu_sample1 <- vector("numeric", length = B)
  mu_sample1[1L] <- rnorm(1L)
  theta_sample1 <- matrix(nrow = B, ncol = 2L)
  theta_sample1[1L, ] <- c(
    runif(1L, min = min(x1), max = max(x1)),
    runif(1L, min = min(x2), max = max(x2))
  )
  acceptace1 <- vector("logical", length = B)
  acceptace1[1L] <- FALSE
  for (i in seq_len(B)[-1L]) {
    # Sample proposal value
    prop <- rmvnorm(1L,
      mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
      sigma = c(10, 0.2, 0.2) * diag(3L)
    )
    mu_p <- prop[1L]
    theta_p <- prop[2L:3L]
    # Compute log ratio of unnormailzed posterior densities
    logr <- log_posterior_unnormalized(mu_p, theta_p, x1, x2) -
      log_posterior_unnormalized(
        mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2
      )
    # Sample uniform random variable
    u <- runif(1L)
    # Accept or reject
    if (isTRUE(log(u) < logr)) {
      mu_sample1[i] <- mu_p
      theta_sample1[i, ] <- theta_p
      acceptace1[i] <- TRUE
    } else {
      mu_sample1[i] <- mu_sample1[i - 1L]
      theta_sample1[i, ] <- theta_sample1[i - 1L, ]
      acceptace1[i] <- FALSE
    }
  }

  # MCMC chain 2
  mu_sample2 <- vector("numeric", length = B)
  mu_sample2[1L] <- rnorm(1L)
  theta_sample2 <- matrix(nrow = B, ncol = 2L)
  theta_sample2[1L, ] <- c(
    runif(1L, min = min(x1), max = max(x1)),
    runif(1L, min = min(x2), max = max(x2))
  )
  acceptace2 <- vector("logical", length = B)
  acceptace2[1L] <- FALSE
  for (i in seq_len(B)[-1L]) {
    # Sample proposal value
    prop <- rmvnorm(1L,
      mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
      sigma = c(10, 0.2, 0.2) * diag(3L)
    )
    mu_p <- prop[1L]
    theta_p <- prop[2L:3L]
    # Compute log ratio of unnormailzed posterior densities
    logr <- log_posterior_unnormalized(mu_p, theta_p, x1, x2) -
      log_posterior_unnormalized(
        mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2
      )
    # Sample uniform random variable
    u <- runif(1L)
    # Accept or reject
    if (isTRUE(log(u) < logr)) {
      mu_sample2[i] <- mu_p
      theta_sample2[i, ] <- theta_p
      acceptace2[i] <- TRUE
    } else {
      mu_sample2[i] <- mu_sample2[i - 1L]
      theta_sample2[i, ] <- theta_sample2[i - 1L, ]
      acceptace2[i] <- FALSE
    }
  }

  # KL divergence
  df <- data.frame(
    index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
    theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
    theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
    acceptace = c(acceptace1, acceptace2),
    chain = rep(c(1L, 2L), each = B)
  )
  pd <- density(df$mu)
  pd_estimate <- approxfun(pd$x, pd$y,
    yleft = .Machine$double.eps,
    yright = .Machine$double.eps
  )
  integrand <- function(x) {
    pd_estimate(x) * {
      log(pd_estimate(x)) - dnorm(x, mean = 0, sd = sd_mu, log = TRUE)
    }
  }
  kl <- tryCatch(integrate(integrand, lower = -Inf, upper = Inf)$value,
    warning = \(x) NA, error = \(x) NA
  )
  # Compute KL by taking the average when integration does not work
  if (is.na(kl)) {
    kl <- mean(
      log(pd_estimate(df$mu)) - dnorm(df$mu, mean = 0, sd = sd_mu, log = TRUE)
    )
  }
  # Potential scale reduction factors
  mu_c1 <- mcmc(df$mu[seq_len(B)])
  mu_c2 <- mcmc(df$mu[seq(B + 1L, 2L * B)])
  theta1_c1 <- mcmc(df$theta1[seq_len(B)])
  theta1_c2 <- mcmc(df$theta1[seq(B + 1L, 2L * B)])
  theta2_c1 <- mcmc(df$theta2[seq_len(B)])
  theta2_c2 <- mcmc(df$theta2[seq(B + 1L, 2L * B)])
  theta1_psrf <- gelman.diag(mcmc.list(theta1_c1, theta1_c2))$psrf[1L]
  mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]
  theta2_psrf <- gelman.diag(mcmc.list(theta2_c1, theta2_c2))$psrf[1L]

  c(kl, mu_psrf, theta1_psrf, theta2_psrf, mean(df$acceptace))
}
stopCluster(cl)


## 6. Results
colnames(result) <-
  c("kl", "mu_psrf", "theta1_psrf", "theta2_psrf", "acceptace")
# Remove invalid values
result <- result[is.finite(result[, "mu_psrf"]), ]
cat("Expected KL:     ",
  mean(result[, "kl"], na.rm = TRUE), " (",
  round(sd(result[, "kl"], na.rm = TRUE) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("mu psrf:         ",
  mean(result[, "mu_psrf"], na.rm = TRUE), " (",
  round(sd(result[, "mu_psrf"], na.rm = TRUE) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("theta1 psrf:     ",
  mean(result[, "theta1_psrf"], na.rm = TRUE), " (",
  round(sd(result[, "theta1_psrf"], na.rm = TRUE) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("theta2 psrf:     ",
  mean(result[, "theta2_psrf"], na.rm = TRUE), " (",
  round(sd(result[, "theta2_psrf"], na.rm = TRUE) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("Acceptance rate: ",
  mean(result[, "acceptace"], na.rm = TRUE), " (",
  round(sd(result[, "acceptace"], na.rm = TRUE) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
