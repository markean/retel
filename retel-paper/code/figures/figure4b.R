## 1. Load packages
options(warn = -1)
options(scipen = 999)
library(coda)
library(ggplot2)
library(melt)
library(mvtnorm)
library(nloptr)
library(retel)


## 2. Constants
# Sample size for each group
n <- 2L
# Standard deviation of prior
sd_mu <- 10
# MCMC parameters
B <- 5000L
# Tau
tau <- 1
# Optimization
ctrl <- el_control(th = 1000, maxit_l = 50L, nthreads = 1L)


## 3. Functions
f <- function(x, par) {
  x - par
}
# Prior and posterior computation
log_prior <- function(mu, theta) {
  dnorm(mu, mean = 0, sd = sd_mu, log = TRUE) -
    sum(log(pi * (1 + (theta - mu)^2L)))
}
loglr_retel_f <- function(theta, x1, x2, tau = 1) {
  retel(f, x1, theta[1L], mean(x1) - theta[1L], 1, tau) +
    retel(f, x2, theta[2L], mean(x2) - theta[2L], 1, tau)
}
log_posterior_unnormalized_retel_f <- function(mu_p, theta, x1, x2, tau = 1) {
  log_prior(mu_p, theta) + loglr_retel_f(theta, x1, x2, tau)
}
loglr_retel_r <- function(theta, x1, x2, tau = 1) {
  retel(f, x1, theta[1L], mean(x1) - theta[1L], 1, tau, "reduced") +
    retel(f, x2, theta[2L], mean(x2) - theta[2L], 1, tau, "reduced")
}
log_posterior_unnormalized_retel_r <- function(mu_p, theta, x1, x2, tau = 1) {
  log_prior(mu_p, theta) + loglr_retel_r(theta, x1, x2, tau)
}
loglr_el <- function(theta, x1, x2) {
  fit1 <- el_mean(x1, par = theta[1L], control = ctrl)
  fit2 <- el_mean(x2, par = theta[2L], control = ctrl)
  lr1 <- if (conv(fit1)) logLR(fit1) else -10000
  lr2 <- if (conv(fit2)) logLR(fit2) else -10000
  lr1 + lr2
}
log_posterior_unnormalized_el <- function(mu_p, theta, x1, x2) {
  log_prior(mu_p, theta) + loglr_el(theta, x1, x2)
}
loglr_etel <- function(theta, x1, x2) {
  if (theta[1L] > max(x1) || theta[1L] < min(x1)) {
    loglr1 <- -10000
  } else {
    loglr1 <- etel(f, x1, theta[1L])
  }
  if (theta[2L] > max(x2) || theta[2L] < min(x2)) {
    loglr2 <- -10000
  } else {
    loglr2 <- etel(f, x2, theta[2L])
  }
  loglr1 + loglr2
}
log_posterior_unnormalized_etel <- function(mu_p, theta, x1, x2) {
  log_prior(mu_p, theta) + loglr_etel(theta, x1, x2)
}


## 4. Simulations
set.seed(25424)
# Generate data
x1 <- rnorm(n, mean = -3, sd = 1)
x2 <- rnorm(n, mean = 3, sd = 1)
# RETEL_f
sd_proposal_mu <- 10
sd_proposal_theta <- 2
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- rnorm(1L)
theta_sample1 <- matrix(nrow = B, ncol = 2L)
theta_sample1[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1L]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_retel_f(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_retel_f(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2, tau
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
theta_sample2[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1L]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_retel_f(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_retel_f(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2, tau
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
df_retel_f <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "retel_f"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_retel_f$mu[seq_len(B)])
mu_c2 <- mcmc(df_retel_f$mu[seq(B + 1L, 2L * B)])
retel_f_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# RETEL_r
sd_proposal_mu <- 10
sd_proposal_theta <- 1
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- rnorm(1L)
theta_sample1 <- matrix(nrow = B, ncol = 2L)
theta_sample1[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace1 <- vector("logical", length = B)
acceptace1[1L] <- FALSE
for (i in seq_len(B)[-1L]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample1[i - 1L], theta_sample1[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_retel_r(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_retel_r(
      mu_sample1[i - 1L], theta_sample1[i - 1L, ], x1, x2, tau
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
theta_sample2[1L, ] <- c(rnorm(1L, mean = -2), rnorm(1L, mean = 2))
acceptace2 <- vector("logical", length = B)
acceptace2[1L] <- FALSE
for (i in seq_len(B)[-1L]) {
  # Sample proposal value
  prop <- rmvnorm(1L,
    mean = c(mu_sample2[i - 1L], theta_sample2[i - 1L, ]),
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_retel_r(mu_p, theta_p, x1, x2, tau) -
    log_posterior_unnormalized_retel_r(
      mu_sample2[i - 1L], theta_sample2[i - 1L, ], x1, x2, tau
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
df_retel_r <- data.frame(
  index = seq_len(2L * B), mu = c(mu_sample1, mu_sample2),
  theta1 = c(theta_sample1[, 1L], theta_sample2[, 1L]),
  theta2 = c(theta_sample1[, 2L], theta_sample2[, 2L]),
  acceptace = c(acceptace1, acceptace2),
  chain = rep(c(1L, 2L), each = B),
  type = "retel_r"
)
# Potential scale reduction factors
mu_c1 <- mcmc(df_retel_r$mu[seq_len(B)])
mu_c2 <- mcmc(df_retel_r$mu[seq(B + 1L, 2L * B)])
retel_r_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# EL
sd_proposal_mu <- 10
sd_proposal_theta <- 0.05
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
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_el(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_el(
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
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_el(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_el(
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
mu_c2 <- mcmc(df_el$mu[seq(B + 1L, 2L * B)])
el_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]

# ETEL
sd_proposal_mu <- 6
sd_proposal_theta <- 0.2
# MCMC chain 1
mu_sample1 <- vector("numeric", length = B)
mu_sample1[1L] <- runif(1L, min = min(x1), max = max(x1))
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
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_etel(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_etel(
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
mu_sample2[1L] <- runif(1L, min = min(x2), max = max(x2))
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
    sigma = c(sd_proposal_mu, sd_proposal_theta, sd_proposal_theta) * diag(3L)
  )
  mu_p <- prop[1L]
  theta_p <- prop[2L:3L]
  # Compute log ratio of unnormailzed posterior densities
  logr <- log_posterior_unnormalized_etel(mu_p, theta_p, x1, x2) -
    log_posterior_unnormalized_etel(
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
mu_c2 <- mcmc(df_etel$mu[seq(B + 1L, 2L * B)])
etel_mu_psrf <- gelman.diag(mcmc.list(mu_c1, mu_c2))$psrf[1L]


# Plot
# 4 x 3
# Potential scale reduction factor
c(retel_f_mu_psrf, retel_r_mu_psrf, el_mu_psrf, etel_mu_psrf)
# Acceptance rate
c(
  mean(df_retel_f$acceptace), mean(df_retel_r$acceptace), mean(df_el$acceptace),
  mean(df_etel$acceptace)
)
legend_labels <- c(
  expression(RETEL[italic(f)]), expression(RETEL[italic(r)]), "EL", "ETEL"
)
df <- rbind(df_retel_f, df_retel_r, df_el, df_etel)
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
    values = c("black", "red", "blue", "green"),
    breaks = c("retel_f", "retel_r", "el", "etel")
  ) +
  scale_linetype_manual(
    labels = legend_labels,
    values = c("solid", "dashed", "dotdash", "twodash"),
    breaks = c("retel_f", "retel_r", "el", "etel")
  ) +
  scale_shape_manual(
    labels = legend_labels,
    values = c(15, 16, 2, 5),
    breaks = c("retel_f", "retel_r", "el", "etel")
  ) +
  xlim(c(-6, 6))
