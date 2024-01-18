# Load packages
library(ggplot2)
library(ggpubr)
library(nloptr)
library(retel)
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-06)

# Data
x <- 0

# Functions
f <- function(x, par) {
  x - par
}

# Grid
grid <- seq(from = -1, to = 1, length.out = 300L)
f_keep_pc <- function(tau) {
  vapply(grid, function(k) {
    retel(f, x,
      par = k, mu = mean(x) - k, Sigma = 1, tau = tau, type = "full",
      opts = opts
    )
  },
  FUN.VALUE = numeric(1L)
  )
}
f_drop_pc <- function(tau) {
  vapply(grid, function(k) {
    retel(f, x,
      par = k, mu = mean(x) - k, Sigma = 1, tau = tau, type = "reduced",
      opts = opts
    )
  },
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
legend_labels <- c(expression(RETEL[italic(f)]), expression(RETEL[italic(r)]))

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
    labels = legend_labels,
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    labels = legend_labels,
    values = c("solid", "dotted"),
    breaks = c("RETEL1", "RETEL2"),
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
    labels = legend_labels,
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    labels = legend_labels,
    values = c("solid", "dotted"),
    breaks = c("RETEL1", "RETEL2")
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
    labels = legend_labels,
    values = c("blue", "red"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_linetype_manual(
    labels = legend_labels,
    values = c("solid", "dotted"),
    breaks = c("RETEL1", "RETEL2")
  ) +
  scale_y_continuous(breaks = c(0, -0.25, -0.5))

# Plots
legend <- get_legend(p1 + theme(
  legend.position = "bottom",
  legend.margin = margin(t = -10),
))
# 8 x 3
ggarrange(p1, p2, p3,
  ncol = 3L, nrow = 1L, common.legend = TRUE, legend = "bottom"
)
