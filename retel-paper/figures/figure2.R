# Load packages
library(ggplot2)
library(ggpubr)
library(nloptr)
library(retel)
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)

# Data
n <- 1
x <- 0

# Functions
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

# Grid
grid <- seq(from = -1, to = 1, length.out = 300)
f_keep_pc <- function(tau) {
  vapply(grid, function(k) RETEL1(x, k, tau, mu = mean(x) - k),
    FUN.VALUE = numeric(1L)
  )
}
f_drop_pc <- function(tau) {
  vapply(grid, function(k) RETEL2(x, k, tau, mu = mean(x) - k),
    FUN.VALUE = numeric(1L)
  )
}
df <- data.frame(
  grid = grid,
  f_keep_pc_1 = f_keep_pc(1),
  f_keep_pc_5 = f_keep_pc(5),
  f_keep_pc_25 = f_keep_pc(25),
  f_drop_pc_1 = f_drop_pc(1),
  f_drop_pc_5 = f_drop_pc(5),
  f_drop_pc_25 = f_drop_pc(25)
)

# tau = 1
p1 <- ggplot(df) +
  geom_line(aes(grid, f_keep_pc_1, color = "RETEL1", linetype = "RETEL1")) +
  geom_line(aes(grid, f_drop_pc_1, color = "RETEL2", linetype = "RETEL2")) +
  geom_vline(xintercept = mean(x), linetype = "dashed") +
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
    x = expression(theta), y = expression(logR),
    title = expression(italic(tau) * " = 1"),
    colour = "", linetype = ""
  ) +
  scale_colour_manual(
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    breaks = c("RETEL1", "RETEL2"),
    values = c("solid", "dotted")
  ) +
  scale_y_continuous(breaks = c(0, -0.25, -0.5))

# tau = 5
p2 <- ggplot(df) +
  geom_line(aes(grid, f_keep_pc_5, color = "RETEL1", linetype = "RETEL1")) +
  geom_line(aes(grid, f_drop_pc_5, color = "RETEL2", linetype = "RETEL2")) +
  geom_vline(xintercept = mean(x), linetype = "dashed") +
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
    x = expression(theta), y = expression(logR),
    title = expression(italic(tau) * " = 5"),
    colour = "", linetype = ""
  ) +
  scale_colour_manual(
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    breaks = c("RETEL1", "RETEL2"),
    values = c("solid", "dotted")
  ) +
  scale_y_continuous(breaks = c(0, -0.25, -0.5))

# tau = 25
p3 <- ggplot(df) +
  geom_line(aes(grid, f_keep_pc_25, color = "RETEL1", linetype = "RETEL1")) +
  geom_line(aes(grid, f_drop_pc_25, color = "RETEL2", linetype = "RETEL2")) +
  geom_vline(xintercept = mean(x), linetype = "dashed") +
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
    x = expression(theta), y = expression(logR),
    title = expression(italic(tau) * " = 25"),
    colour = "", linetype = ""
  ) +
  scale_colour_manual(
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    breaks = c("RETEL1", "RETEL2"),
    values = c("solid", "dotted")
  ) +
  scale_y_continuous(breaks = c(0, -0.25, -0.5))

# Plots
legend <- get_legend(p1 + theme(
  legend.position = "bottom",
  legend.margin = margin(t = -10),
))
ggarrange(p1, p2, p3,
  ncol = 3L, nrow = 1L, common.legend = TRUE,
  legend = "bottom"
)
