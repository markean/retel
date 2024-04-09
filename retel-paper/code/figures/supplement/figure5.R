## 1. Load packages
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(melt)
library(mvtnorm)
library(nloptr)
library(retel)


## 2. Constants
p <- 10
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)
S <- 1e+03L


## 3. Functions
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


## 4. Simulation
# n = 20 ----
n <- 20
an <- log(n) / 2
par <- rep(0, p)
set.seed(6325203)
n20_res_retel_f <- vector("numeric", length = S)
n20_res_retel_r <- vector("numeric", length = S)
n20_res_aetel <- vector("numeric", length = S)
n20_res_ael <- vector("numeric", length = S)
n20_res_bael <- vector("numeric", length = S)
for (i in seq_len(S)) {
  # Data
  x <- rmvnorm(n, sigma = diag(p))

  # RETEL_f
  n20_res_retel_f[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "full"
  )

  # RETEL_r
  n20_res_retel_r[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "reduced"
  )

  # AETEL
  n20_res_aetel[i] <- -2 * aetel(x, par)

  # AEL
  p_ael <- par - an * (colMeans(x) - par)
  n20_res_ael[i] <- el_mean(rbind(x, p_ael), par = par) |>
    chisq()

  # BAEL
  v <- colMeans(x) - par
  u <- v / norm(v, type = "2")
  p1_bael <- par - 1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  p2_bael <- 2 * colMeans(x) - par +
    1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  n20_res_bael[i] <- el_mean(rbind(x, p1_bael, p2_bael), par = par) |>
    chisq()
}

# n = 50 ----
n <- 50
an <- log(n) / 2
par <- rep(0, p)
cv <- qchisq(0.95, df = p)
set.seed(6325203)
n50_res_retel_f <- vector("numeric", length = S)
n50_res_retel_r <- vector("numeric", length = S)
n50_res_aetel <- vector("numeric", length = S)
n50_res_ael <- vector("numeric", length = S)
n50_res_bael <- vector("numeric", length = S)
for (i in seq_len(S)) {
  # Data
  x <- rmvnorm(n, sigma = diag(p))

  # RETEL_f
  n50_res_retel_f[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "full"
  )

  # RETEL_r
  n50_res_retel_r[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "reduced"
  )

  # AETEL
  n50_res_aetel[i] <- -2 * aetel(x, par)

  # AEL
  p_ael <- par - an * (colMeans(x) - par)
  n50_res_ael[i] <- el_mean(rbind(x, p_ael), par = par) |>
    chisq()

  # BAEL
  v <- colMeans(x) - par
  u <- v / norm(v, type = "2")
  p1_bael <- par - 1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  p2_bael <- 2 * colMeans(x) - par +
    1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  n50_res_bael[i] <- el_mean(rbind(x, p1_bael, p2_bael), par = par) |>
    chisq()
}

# n = 100 ----
n <- 100
an <- log(n) / 2
par <- rep(0, p)
cv <- qchisq(0.95, df = p)
set.seed(6325203)
n100_res_retel_f <- vector("numeric", length = S)
n100_res_retel_r <- vector("numeric", length = S)
n100_res_aetel <- vector("numeric", length = S)
n100_res_ael <- vector("numeric", length = S)
n100_res_bael <- vector("numeric", length = S)
for (i in seq_len(S)) {
  # Data
  x <- rmvnorm(n, sigma = diag(p))

  # RETEL_f
  n100_res_retel_f[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "full"
  )

  # RETEL_r
  n100_res_retel_r[i] <- -2 * retel(f, x, par,
    mu = colMeans(x),
    Sigma = sqrt(p) * cov(x),
    tau = log(n),
    type = "reduced"
  )

  # AETEL
  n100_res_aetel[i] <- -2 * aetel(x, par)

  # AEL
  p_ael <- par - an * (colMeans(x) - par)
  n100_res_ael[i] <- el_mean(rbind(x, p_ael), par = par) |>
    chisq()

  # BAEL
  v <- colMeans(x) - par
  u <- v / norm(v, type = "2")
  p1_bael <- par - 1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  p2_bael <- 2 * colMeans(x) - par +
    1.9 * u / sqrt(colSums(u * (solve(cov(x)) %*% u)))
  n100_res_bael[i] <- el_mean(rbind(x, p1_bael, p2_bael), par = par) |>
    chisq()
}


## 5. Plots
# n = 20 ----
p_n20_retel_f <- ggplot(
  data.frame(sample = n20_res_retel_f),
  aes(sample = n20_res_retel_f)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  )

p_n20_retel_r <- ggplot(
  data.frame(sample = n20_res_retel_r),
  aes(sample = n20_res_retel_r)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  )

p_n20_aetel <- ggplot(
  data.frame(sample = n20_res_aetel),
  aes(sample = n20_res_aetel)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "",
    linetype = ""
  )

p_n20_ael <- ggplot(
  data.frame(sample = n20_res_ael),
  aes(sample = n20_res_ael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AEL), colour = "",
    linetype = ""
  )

p_n20_bael <- ggplot(
  data.frame(sample = n20_res_bael),
  aes(sample = n20_res_bael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(BAEL), colour = "",
    linetype = ""
  )

# n = 50 ----
p_n50_retel_f <- ggplot(
  data.frame(sample = n50_res_retel_f),
  aes(sample = n50_res_retel_f)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  )

p_n50_retel_r <- ggplot(
  data.frame(sample = n50_res_retel_r),
  aes(sample = n50_res_retel_r)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  )

p_n50_aetel <- ggplot(
  data.frame(sample = n50_res_aetel),
  aes(sample = n50_res_aetel)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "",
    linetype = ""
  )

p_n50_ael <- ggplot(
  data.frame(sample = n50_res_ael),
  aes(sample = n50_res_ael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AEL), colour = "",
    linetype = ""
  )

p_n50_bael <- ggplot(
  data.frame(sample = n50_res_bael),
  aes(sample = n50_res_bael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(BAEL), colour = "",
    linetype = ""
  )

# n = 100 ----
p_n100_retel_f <- ggplot(
  data.frame(sample = n100_res_retel_f),
  aes(sample = n100_res_retel_f)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  )

p_n100_retel_r <- ggplot(
  data.frame(sample = n100_res_retel_r),
  aes(sample = n100_res_retel_r)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  )

p_n100_aetel <- ggplot(
  data.frame(sample = n100_res_aetel),
  aes(sample = n100_res_aetel)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "",
    linetype = ""
  )

p_n100_ael <- ggplot(
  data.frame(sample = n100_res_ael),
  aes(sample = n100_res_ael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AEL), colour = "",
    linetype = ""
  )

p_n100_bael <- ggplot(
  data.frame(sample = n100_res_bael),
  aes(sample = n100_res_bael)
) +
  stat_qq(distribution = qchisq, dparams = list(df = p), size = 0.5) +
  stat_qq_line(distribution = qchisq, dparams = list(df = p)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(BAEL), colour = "",
    linetype = ""
  )

plots <- plot_grid(
  p_n20_retel_f, p_n20_retel_r, p_n20_aetel, p_n20_ael, p_n20_bael,
  p_n50_retel_f, p_n50_retel_r, p_n50_aetel, p_n50_ael, p_n50_bael,
  p_n100_retel_f, p_n100_retel_r, p_n100_aetel, p_n100_ael, p_n100_bael,
  nrow = 3L, align = "v"
)
y.grob <- textGrob("Empirical quantiles", gp = gpar(fontsize = 12), rot = 90)
x.grob <- textGrob(
  expression(italic(chi^2) * italic((p)) * " quantiles"),
  gp = gpar(fontsize = 12)
)
# 850 x 700 (PNG)
grid.arrange(arrangeGrob(plots, left = y.grob, bottom = x.grob))
