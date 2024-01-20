# Load packages
library(ggplot2)
suppressMessages(library(here))
suppressMessages(i_am("code/figures/figure4a.R"))

# Data
df <- data.frame(
  method = rep(c("retel_f", "retel_r", "el", "etel"), each = 5L),
  n = rep(c(2L, 4L, 6L, 8L, 10L), 4L)
)
ekl_retel_f <- c(
  readRDS(here("simulations/kl/retel_f/n2.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n4.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n6.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n8.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n10.rds"))[, "kl"] |>
    mean(na.rm = TRUE)
)
ekl_retel_r <- c(
  readRDS(here("simulations/kl/retel_r/n2.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n4.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n6.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n8.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n10.rds"))[, "kl"] |>
    mean(na.rm = TRUE)
)
ekl_el <- c(
  readRDS(here("simulations/kl/el/n2.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n4.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n6.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n8.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n10.rds"))[, "kl"] |>
    mean(na.rm = TRUE)
)
ekl_etel <- c(
  readRDS(here("simulations/kl/etel/n2.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n4.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n6.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n8.rds"))[, "kl"] |>
    mean(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n10.rds"))[, "kl"] |>
    mean(na.rm = TRUE)
)
se_retel_f <- (c(
  readRDS(here("simulations/kl/retel_f/n2.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n4.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n6.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n8.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_f/n10.rds"))[, "kl"] |>
    sd(na.rm = TRUE)
) / sqrt(1000L)) |>
  round(4L)
se_retel_r <- (c(
  readRDS(here("simulations/kl/retel_r/n2.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n4.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n6.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n8.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/retel_r/n10.rds"))[, "kl"] |>
    sd(na.rm = TRUE)
) / sqrt(1000L)) |>
  round(4L)
se_el <- (c(
  readRDS(here("simulations/kl/el/n2.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n4.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n6.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n8.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/el/n10.rds"))[, "kl"] |>
    sd(na.rm = TRUE)
) / sqrt(1000L)) |>
  round(4L)
se_etel <- (c(
  readRDS(here("simulations/kl/etel/n2.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n4.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n6.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n8.rds"))[, "kl"] |>
    sd(na.rm = TRUE),
  readRDS(here("simulations/kl/etel/n10.rds"))[, "kl"] |>
    sd(na.rm = TRUE)
) / sqrt(1000L)) |>
  round(4L)
df$ekl <- c(ekl_retel_f, ekl_retel_r, ekl_el, ekl_etel)
df$se <- c(se_retel_f, se_retel_r, se_el, se_etel)

# Plot
# 4 x 3
legend_labels <- c(
  expression(RETEL[italic(f)]), expression(RETEL[italic(r)]), "EL", "ETEL"
)
ggplot(df, aes(n, ekl, color = method, group = method, linetype = method)) +
  geom_line() +
  geom_point(aes(shape = method)) +
  geom_errorbar(aes(ymin = ekl - se, ymax = ekl + se),
    width = 0.2,
    position = position_dodge(0.15)
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
    x = expression(n), y = expression(E * "[KL]"),
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
  )
