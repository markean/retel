# Load packages
library(ggplot2)

# Data
df <- data.frame(
  method = rep(c("retel1", "retel2", "el", "etel"), each = 5L),
  n = rep(c(2L, 4L, 6L, 8L, 10L), 4L)
)
kl_retel1 <- c(0.8609873, 0.8721041, 0.8856395, 0.8961443, 0.8989103)
kl_retel2 <- c(0.8636815, 0.8822516, 0.8865799, 0.8969079, 0.9007792)
kl_el <- c(0.9144397, 0.8846304, 0.8832214, 0.8858823, 0.8842803)
kl_etel <- c(0.9173491, 0.8884114, 0.8868945, 0.8827214, 0.8833329)
se_retel1 <- c(0.0039, 0.0024, 0.0023, 0.0021, 0.0021)
se_retel2 <- c(0.0038, 0.0025, 0.0022, 0.0021, 0.0022)
se_el <- c(0.0045, 0.0025, 0.0021, 0.0019, 0.0017)
se_etel <- c(0.0056, 0.0026, 0.0022, 0.0019, 0.0016)
df$kl <- c(kl_retel1, kl_retel2, kl_el, kl_etel)
df$se <- c(se_retel1, se_retel2, se_el, se_etel)

# Plot
legend_labels <- c(
  expression(RETEL[italic(f)]), expression(RETEL[italic(r)]), "EL", "ETEL"
)
ggplot(df, aes(n, kl, color = method, group = method, linetype = method)) +
  geom_line() +
  geom_point(aes(shape = method)) +
  geom_errorbar(aes(ymin = kl - se, ymax = kl + se),
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
  )
