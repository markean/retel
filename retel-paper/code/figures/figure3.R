# Load packages
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
suppressMessages(library(here))
suppressMessages(here::i_am("code/figures/figure3.R"))

# Data
result <- readRDS(here::here("simulations/mb/mb_1/n5s5.rds"))
df <- as.data.frame(result)

# RETEL_f
p1 <- ggplot(df, aes(sample = retel_f)) +
  stat_qq(distribution = qunif, size = 0.5) +
  stat_qq_line(distribution = qunif) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    # Reduce the space between x axis label and legend
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]),
    # Remove excessive legend
    colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# RETEL_r
p2 <- ggplot(df, aes(sample = retel_r)) +
  stat_qq(distribution = qunif, size = 0.5) +
  stat_qq_line(distribution = qunif) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    # Reduce the space between x axis label and legend
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]),
    # Remove excessive legend
    colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# ETEL
p3 <- ggplot(df, aes(sample = etel)) +
  stat_qq(distribution = qunif, size = 0.5) +
  stat_qq_line(distribution = qunif) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    # Reduce the space between x axis label and legend
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(ETEL),
    # Remove excessive legend
    colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# AETEL
p4 <- ggplot(df, aes(sample = aetel)) +
  stat_qq(distribution = qunif, size = 0.5) +
  stat_qq_line(distribution = qunif) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = alpha("white", 0)),
    # Reduce the space between x axis label and legend
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL),
    # Remove excessive legend
    colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# Plots
plots <- plot_grid(p1, p2, p3, p4, nrow = 1L, align = "v")
y.grob <- textGrob(expression(italic(H) * " quantiles"),
  gp = gpar(fontsize = 12), rot = 90
)
x.grob <- textGrob(expression(italic(U(0, 1)) * " quantiles"),
  gp = gpar(fontsize = 12)
)
# 8 x 3
grid.arrange(arrangeGrob(plots, left = y.grob, bottom = x.grob))
