# Load packages
library(ggplot2)
library(nloptr)
library(retel)
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)

# Data
x <- c(-2, 2)

# Functions
f <- function(x, par) {
  x - par
}
discrete_obj <- function(l, par, pg) {
  sum(exp(l * (x - par))) + mean(exp(l * pg))
}
gr_discrete_obj <- function(l, par, pg) {
  sum(exp(l * (x - par)) * (x - par)) + mean(exp(l * pg) * pg)
}
lambda_WETEL <- function(par, m) {
  pg <- qnorm(1:m / (m + 1), mean = 0, sd = 1)
  out <- nloptr(
    x0 = 0, eval_f = discrete_obj, eval_grad_f = gr_discrete_obj, opts = opts,
    par = par, pg = pg
  )
  out$solution
}
lambda_RETEL <- function(par) {
  out <- retel(f, x, par, mu = 0, Sigma = 1, tau = 1, opts = opts)
  attr(out, "optim")$solution
}

# Convex hull constraint satisfied at 1
target <- lambda_RETEL(1)
grid <- seq(from = 5, to = 20, length.out = 100)
approx <- vapply(round(2^grid), function(k) lambda_WETEL(1, k),
  FUN.VALUE = numeric(1L)
)
df <- data.frame(m = grid, lambda = approx)
# 4 x 3
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
# 4 x 3
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
