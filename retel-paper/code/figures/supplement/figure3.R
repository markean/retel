## Load packages
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)


## Load data
df_taulogn_n5s1 <-
  readRDS("./simulations/mb/mb_logn/n5s1.rds") |>
  as.data.frame()
df_taulogn_n20s1 <-
  readRDS("./simulations/mb/mb_logn/n20s1.rds") |>
  as.data.frame()
df_taulogn_n50s1 <-
  readRDS("./simulations/mb/mb_logn/n50s1.rds") |>
  as.data.frame()
df_taulogn_n100s1 <-
  readRDS("./simulations/mb/mb_logn/n100s1.rds") |>
  as.data.frame()


## n = 5 ----
# RETEL_f
p_taulogn_n5s1_1 <- ggplot(df_taulogn_n5s1, aes(sample = retel_f)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# RETEL_r
p_taulogn_n5s1_2 <- ggplot(df_taulogn_n5s1, aes(sample = retel_r)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# ETEL
p_taulogn_n5s1_3 <- ggplot(df_taulogn_n5s1, aes(sample = etel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(ETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# AETEL
p_taulogn_n5s1_4 <- ggplot(df_taulogn_n5s1, aes(sample = aetel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


## n = 20 ----
# RETEL_f
p_taulogn_n20s1_1 <- ggplot(df_taulogn_n20s1, aes(sample = retel_f)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# RETEL_r
p_taulogn_n20s1_2 <- ggplot(df_taulogn_n20s1, aes(sample = retel_r)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# ETEL
p_taulogn_n20s1_3 <- ggplot(df_taulogn_n20s1, aes(sample = etel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(ETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# AETEL
p_taulogn_n20s1_4 <- ggplot(df_taulogn_n20s1, aes(sample = aetel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


## n = 50 ----
# RETEL_f
p_taulogn_n50s1_1 <- ggplot(df_taulogn_n50s1, aes(sample = retel_f)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# RETEL_r
p_taulogn_n50s1_2 <- ggplot(df_taulogn_n50s1, aes(sample = retel_r)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# ETEL
p_taulogn_n50s1_3 <- ggplot(df_taulogn_n50s1, aes(sample = etel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(ETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# AETEL
p_taulogn_n50s1_4 <- ggplot(df_taulogn_n50s1, aes(sample = aetel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


## n = 100 ----
# RETEL_f
p_taulogn_n100s1_1 <- ggplot(df_taulogn_n100s1, aes(sample = retel_f)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(f)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# RETEL_r
p_taulogn_n100s1_2 <- ggplot(df_taulogn_n100s1, aes(sample = retel_r)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(RETEL[italic(r)]), colour = "",
    linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# ETEL
p_taulogn_n100s1_3 <- ggplot(df_taulogn_n100s1, aes(sample = etel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(ETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))

# AETEL
p_taulogn_n100s1_4 <- ggplot(df_taulogn_n100s1, aes(sample = aetel)) +
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
    legend.margin = margin(t = -10),
    legend.key = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, y = NULL, title = expression(AETEL), colour = "", linetype = ""
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


## Plot ----
plots <- plot_grid(
  p_taulogn_n5s1_1, p_taulogn_n5s1_2, p_taulogn_n5s1_3, p_taulogn_n5s1_4,
  p_taulogn_n20s1_1, p_taulogn_n20s1_2, p_taulogn_n20s1_3, p_taulogn_n20s1_4,
  p_taulogn_n50s1_1, p_taulogn_n50s1_2, p_taulogn_n50s1_3, p_taulogn_n50s1_4,
  p_taulogn_n100s1_1, p_taulogn_n100s1_2, p_taulogn_n100s1_3,
  p_taulogn_n100s1_4,
  nrow = 4L, align = "v"
)
y.grob <- textGrob(expression(italic(H) * " quantiles"),
  gp = gpar(fontsize = 12), rot = 90
)
x.grob <- textGrob(expression(italic(U(0, 1)) * " quantiles"),
  gp = gpar(fontsize = 12)
)
# 850 x 1000 (PNG)
grid.arrange(arrangeGrob(plots, left = y.grob, bottom = x.grob))
