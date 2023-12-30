## 1. Load packages ----
library(coda)
library(ggplot2)
library(melt)
library(mvtnorm)
library(nloptr)


## 2. Functions ----
d <- function(l, x, theta) {
  mean(exp(l * (x - theta)))
}
# ETEL
obj <- function(l, x, theta) {
  d(l, x, theta)
}
gr_obj <- function(l, x, theta) {
  mean(exp(l * (x - theta)) * (x - theta))
}
ETEL <- function(x, theta) {
  if (theta > max(x) || theta < min(x)) {
    return(-10000)
  }
  n <- length(x)
  out <- nloptr(
    x0 = 0, eval_f = obj, eval_grad_f = gr_obj, opts = opts, theta = theta,
    x = x,
  )
  lambda <- out$solution
  -n * log(d(lambda, x, theta)) + lambda * sum(x - theta)
}
# RETEL
d2 <- function(x, l, theta, tau) {
  n <- length(x)
  sum(exp(l * (x - theta))) / (n + tau)
}
penalty2 <- function(l, tau, mu = 0, sigma = 1) {
  (tau / (n + tau)) * exp(l * mu + 0.5 * l^2 * sigma^2)
}
continuous_obj2 <- function(l, x, theta, tau, mu = 0, sigma = 1) {
  d2(x, l, theta, tau) + penalty2(l, tau, mu, sigma)
}
gr_continuous_obj2 <- function(l, x, theta, tau, mu = 0, sigma = 1) {
  n <- length(x)
  sum(exp(l * (x - theta)) * (x - theta)) / (n + tau) +
    (tau / (n + tau)) * exp(l * mu + 0.5 * l^2 * sigma^2) * (mu + sigma^2 * l)
}
RETEL1 <- function(x, theta, tau = 1, mu = 0, sigma = 1) {
  n <- length(x)
  out <- nloptr(
    x0 = 0, eval_f = continuous_obj2,
    eval_grad_f = gr_continuous_obj2, opts = opts, theta = theta, tau = tau,
    x = x, mu = mu, sigma = sigma
  )
  lambda <- out$solution
  -(n + 1) * log(continuous_obj2(lambda, x, theta, tau, mu, sigma)) +
    lambda * sum(x - theta) + (lambda * mu + 0.5 * lambda^2 * sigma^2)
}
RETEL2 <- function(x, theta, tau = 1, mu = 0, sigma = 1) {
  n <- length(x)
  out <- nloptr(
    x0 = 0, eval_f = continuous_obj2,
    eval_grad_f = gr_continuous_obj2, opts = opts, theta = theta, tau = tau,
    x = x, mu = mu, sigma = sigma
  )
  l <- out$solution
  -n * log(continuous_obj2(l, x, theta, tau, mu, sigma)) + l * sum(x - theta)
}
# Prior and posterior computation
log_prior <- function(mu, theta) {
  dnorm(mu, mean = 0, sd = sd_mu, log = TRUE) -
    sum(log(pi * (1 + (theta - mu)^2)))
}
logLR_EL <- function(theta, x1, x2) {
  fit1 <- el_mean(x1, par = theta[1L], control = ctrl)
  fit2 <- el_mean(x2, par = theta[2L], control = ctrl)
  lr1 <- if (conv(fit1)) logLR(fit1) else -10000
  lr2 <- if (conv(fit2)) logLR(fit2) else -10000
  lr1 + lr2
}
log_posterior_unnormalized_EL <- function(mu_p, theta, x1, x2) {
  log_prior(mu_p, theta) + logLR_EL(theta, x1, x2)
}

logLR_ETEL <- function(theta, x1, x2) {
  ETEL(x1, theta[1L]) + ETEL(x2, theta[2L])
}
log_posterior_unnormalized_ETEL <- function(mu_p, theta, x1, x2) {
  log_prior(mu_p, theta) + logLR_ETEL(theta, x1, x2)
}
logLR_RETEL1 <- function(theta, x1, x2, tau = 1) {
  RETEL1(x1, theta[1L], tau, mu = mean(x1) - theta[1L], sigma = 1) +
    RETEL1(x2, theta[2L], tau, mu = mean(x2) - theta[2L], sigma = 1)
}
log_posterior_unnormalized_RETEL1 <- function(mu_p, theta, x1, x2, tau = 1) {
  log_prior(mu_p, theta) + logLR_RETEL1(theta, x1, x2, tau)
}
logLR_RETEL2 <- function(theta, x1, x2, tau = 1) {
  RETEL2(x1, theta[1L], tau, mu = mean(x1) - theta[1L], sigma = 1) +
    RETEL2(x2, theta[2L], tau, mu = mean(x2) - theta[2L], sigma = 1)
}
log_posterior_unnormalized_RETEL2 <- function(mu_p, theta, x1, x2, tau = 1) {
  log_prior(mu_p, theta) + logLR_RETEL2(theta, x1, x2, tau)
}


## 3. Parameters ----
# Sample size for each group
n <- 2
# Standard deviation of prior
sd_mu <- 10
# MCMC parameters
B <- 5000
# Tau
tau <- 1
# Optimization
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
ctrl <- el_control(th = 1000, maxit_l = 50, nthreads = 1L)


## 4. Simulation ----
set.seed(25424)
# Generate data
x1 <- rnorm(n, mean = -3, sd = 1)
x2 <- rnorm(n, mean = 3, sd = 1)
# RETEL_f ----
sd_proposal_mu <- 10
sd_proposal_theta <- 2
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- rnorm(1)
theta_sample1 <- matrix(nrow = B, ncol = 2)
theta_sample1[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_RETEL1(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_RETEL1(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2, tau
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample1[i] <- mu_p
    theta_sample1[i, ] <- theta_p
    acceptace1[i] <- TRUE
  } else {
    mu_sample1[i] <- mu_sample1[i - 1L]
    theta_sample1[i, ] <- theta_sample1[i - 1, ]
    acceptace1[i] <- FALSE
  }
}
# MCMC chain 2
mu_sample2 <- vector("numeric", length = B)
mu_sample2[1L] <- rnorm(1)
theta_sample2 <- matrix(nrow = B, ncol = 2)
theta_sample2[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_RETEL1(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_RETEL1(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2, tau
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample2[i] <- mu_p
    theta_sample2[i, ] <- theta_p
    acceptace2[i] <- TRUE
  } else {
    mu_sample2[i] <- mu_sample2[i - 1L]
    theta_sample2[i, ] <- theta_sample2[i - 1, ]
    acceptace2[i] <- FALSE
  }
}
# KL divergence
df_retel1 <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "retel1"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_retel1$mu[seq_len(B)])
mu_c2 <- mcmc(df_retel1$mu[seq(B + 1, 2 * B)])
retel1_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# RETEL_r ----
sd_proposal_mu <- 10
sd_proposal_theta <- 1
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- rnorm(1)
theta_sample1 <- matrix(nrow = B, ncol = 2)
theta_sample1[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_RETEL2(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_RETEL2(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2, tau
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample1[i] <- mu_p
    theta_sample1[i, ] <- theta_p
    acceptace1[i] <- TRUE
  } else {
    mu_sample1[i] <- mu_sample1[i - 1L]
    theta_sample1[i, ] <- theta_sample1[i - 1, ]
    acceptace1[i] <- FALSE
  }
}
# MCMC chain 2
mu_sample2 <- vector("numeric", length = B)
mu_sample2[1L] <- rnorm(1)
theta_sample2 <- matrix(nrow = B, ncol = 2)
theta_sample2[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_RETEL2(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_RETEL2(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2, tau
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample2[i] <- mu_p
    theta_sample2[i, ] <- theta_p
    acceptace2[i] <- TRUE
  } else {
    mu_sample2[i] <- mu_sample2[i - 1L]
    theta_sample2[i, ] <- theta_sample2[i - 1, ]
    acceptace2[i] <- FALSE
  }
}
# KL divergence
df_retel2 <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "retel2"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_retel2$mu[seq_len(B)])
mu_c2 <- mcmc(df_retel2$mu[seq(B + 1, 2 * B)])
retel2_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# EL ----
sd_proposal_mu <- 10
sd_proposal_theta <- 0.05
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- rnorm(1)
theta_sample1 <- matrix(nrow = B, ncol = 2)
theta_sample1[1L, ] <- c(
  runif(1L, min = min(x1), max = max(x1)),
  runif(1L, min = min(x2), max = max(x2))
)
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_EL(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_EL(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample1[i] <- mu_p
    theta_sample1[i, ] <- theta_p
    acceptace1[i] <- TRUE
  } else {
    mu_sample1[i] <- mu_sample1[i - 1L]
    theta_sample1[i, ] <- theta_sample1[i - 1, ]
    acceptace1[i] <- FALSE
  }
}
# MCMC chain 2
mu_sample2 <- vector("numeric", length = B)
mu_sample2[1L] <- rnorm(1)
theta_sample2 <- matrix(nrow = B, ncol = 2)
theta_sample2[1L, ] <- c(
  runif(1L, min = min(x1), max = max(x1)),
  runif(1L, min = min(x2), max = max(x2))
)
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_EL(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_EL(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample2[i] <- mu_p
    theta_sample2[i, ] <- theta_p
    acceptace2[i] <- TRUE
  } else {
    mu_sample2[i] <- mu_sample2[i - 1L]
    theta_sample2[i, ] <- theta_sample2[i - 1, ]
    acceptace2[i] <- FALSE
  }
}
# KL divergence
df_el <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "el"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_el$mu[seq_len(B)])
mu_c2 <- mcmc(df_el$mu[seq(B + 1, 2 * B)])
el_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# ETEL ----
sd_proposal_mu <- 6
sd_proposal_theta <- 0.2
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- runif(1, min = min(x1), max = max(x1))
theta_sample1 <- matrix(nrow = B, ncol = 2)
theta_sample1[1L, ] <- c(
  runif(1L, min = min(x1), max = max(x1)),
  runif(1L, min = min(x2), max = max(x2))
)
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_ETEL(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_ETEL(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample1[i] <- mu_p
    theta_sample1[i, ] <- theta_p
    acceptace1[i] <- TRUE
  } else {
    mu_sample1[i] <- mu_sample1[i - 1L]
    theta_sample1[i, ] <- theta_sample1[i - 1, ]
    acceptace1[i] <- FALSE
  }
}
# MCMC chain 2
mu_sample2 <- vector("numeric", length = B)
mu_sample2[1L] <- runif(1, min = min(x2), max = max(x2))
theta_sample2 <- matrix(nrow = B, ncol = 2)
theta_sample2[1L, ] <- c(
  runif(1L, min = min(x1), max = max(x1)),
  runif(1L, min = min(x2), max = max(x2))
)
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_ETEL(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_ETEL(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2
    )
  # Sample uniform random variable
  u <- runif(1)
  # Accept or reject
  if (isTRUE(log(u) < logr)) {
    mu_sample2[i] <- mu_p
    theta_sample2[i, ] <- theta_p
    acceptace2[i] <- TRUE
  } else {
    mu_sample2[i] <- mu_sample2[i - 1L]
    theta_sample2[i, ] <- theta_sample2[i - 1, ]
    acceptace2[i] <- FALSE
  }
}
# KL divergence
df_etel <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "etel"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_etel$mu[seq_len(B)])
mu_c2 <- mcmc(df_etel$mu[seq(B + 1, 2 * B)])
etel_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]


# Plot ----
# Potential scale reduction factor
c(retel1_mu_psrf, retel2_mu_psrf, el_mu_psrf, etel_mu_psrf)
# Acceptance rate
c(
  mean(df_retel1$acceptace), mean(df_retel2$acceptace), mean(df_el$acceptace),
  mean(df_etel$acceptace)
)
legend_labels <- c(
  expression(RETEL[italic(f)]), expression(RETEL[italic(r)]), "EL", "ETEL"
)
df <- rbind(df_retel1, df_retel2, df_el, df_etel)
ggplot(df, aes(x = mu, color = type, linetype = type)) +
  geom_density(adjust = 2) +
  geom_vline(
    xintercept = c(min(x1), max(x1)), alpha = 1, linetype = "dotted"
  ) +
  geom_vline(
    xintercept = c(min(x2), max(x2)), alpha = 1, linetype = "dotted"
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = expression(italic(mu)),
    y = expression(pi * "(" * mu * "|" * D[n] * ")"),
    colour = "", linetype = "", shape = ""
  ) +
  scale_colour_manual(
    labels = legend_labels,
    breaks = c("retel1", "retel2", "el", "etel"),
    values = c("black", "red", "blue", "green")
  ) +
  scale_linetype_manual(
    labels = legend_labels,
    breaks = c("retel1", "retel2", "el", "etel"),
    values = c("solid", "dashed", "dotdash", "twodash")
  ) +
  scale_shape_manual(
    labels = legend_labels,
    breaks = c("retel1", "retel2", "el", "etel"),
    values = c(15, 16, 2, 5)
  ) +
  xlim(c(-6, 6))
