# Load packages
library(ggplot2)
library(melt)
library(nloptr)
library(retel)
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)

# Data
x <- c(-2, 2)

# Functions
d <- function(l, theta) {
  sum(exp(l * (x - theta)))
}
penalty_m <- function(l, pg) {
  mean(exp(l * pg))
}
penalty <- function(l) {
  exp(0.5 * l^2)
}
discrete_obj <- function(l, theta, pg) {
  d(l, theta) + penalty_m(l, pg)
}
continuous_obj <- function(l, theta) {
  d(l, theta) + penalty(l)
}
gr_discrete_obj <- function(l, theta, pg) {
  sum(exp(l * (x - theta)) * (x - theta)) + mean(exp(l * pg) * pg)
}
gr_continuous_obj <- function(l, theta) {
  sum(exp(l * (x - theta)) * (x - theta)) + exp(0.5 * l^2) * l
}
lambda_WETEL <- function(theta, m) {
  pg <- qnorm(1:m / (m + 1), mean = 0, sd = 1)
  out <- nloptr(
    x0 = 0, eval_f = discrete_obj, eval_grad_f = gr_discrete_obj, opts = opts,
    theta = theta, pg = pg
  )
  out$solution
}
lambda_RETEL <- function(theta) {
  out <- nloptr(
    x0 = 0, eval_f = continuous_obj, eval_grad_f = gr_continuous_obj,
    opts = opts, theta = theta
  )
  out$solution
}

# Convex hull constraint satisfied at 1
target <- lambda_RETEL(1)
grid <- seq(from = 5, to = 20, length.out = 100)
approx <- vapply(round(2^grid), function(k) lambda_WETEL(1, k),
  FUN.VALUE = numeric(1L)
)
df <- data.frame(m = grid, lambda = approx)
ggplot(df) +
  geom_line(aes(m, lambda)) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  ) +
  labs(
    x = expression(italic(log[2] * "m")), y = expression(italic(lambda[WET])),
  ) +
  geom_hline(yintercept = target, linetype = "dashed") +
  scale_y_continuous(
    breaks = c(target, 0.233, 0.235, 0.237),
    labels = c(
      expression(italic(lambda[RET])),
      "0.233", "0.235", "0.237"
    )
  )

# Convex hull constraint violated at 3
target2 <- lambda_RETEL(-3)
grid <- seq(from = 5, to = 20, length.out = 100)
approx2 <- vapply(round(2^grid), function(k) lambda_WETEL(-3, k),
  FUN.VALUE = numeric(1L)
)
df2 <- data.frame(m = grid, lambda = approx2)
ggplot(df2) +
  geom_line(aes(m, lambda)) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  ) +
  labs(
    x = expression(italic(log[2] * "m")), y = expression(italic(lambda[WET])),
  ) +
  geom_hline(yintercept = target2, linetype = "dashed") +
  scale_y_continuous(
    breaks = c(-0.68, -0.66, -0.64, target2),
    labels = c("-0.68", "-0.66", "-0.64", expression(italic(lambda[RET])))
  )
