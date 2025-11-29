library(ggplot2)
library(tidyr)
library(dplyr)

# Number of states
N <- 1:6

# Black–Scholes HMM parameter counts:
# 2 state-dependent: N^2 - N + 2N = N^2 + N
# 1 state-dependent: N^2 - N + (N + 1) = N^2 + 1
# 0 state-dependent: N^2 - N + 2
data <- data.frame(
  N = N,
  Parameters_2 = N^2 + N,          # both μ and σ state-dependent
  Parameters_1 = N^2 + 1,          # exactly one of μ or σ state-dependent
  Parameters_0 = N^2 - N + 2       # neither state-dependent
)

# Long format
long_data <- pivot_longer(
  data,
  cols = starts_with("Parameters"),
  names_to = "Case",
  values_to = "Parameters"
)

# Relabel legend entries
long_data$Case <- recode(
  long_data$Case,
  Parameters_2 = "2 state-dependent",
  Parameters_1 = "1 state-dependent",
  Parameters_0 = "0 state-dependent"
)

# Shapes per case
shape_map <- c(
  "2 state-dependent" = 15,  # square
  "1 state-dependent"        = 17,  # triangle
  "0 state-dependent"        = 16   # circle
)

# Colors per case
color_map <- c(
  "2 state-dependent" = "#CC5500",
  "1 state-dependent"        = "#404080",
  "0 state-dependent"        = "#87A96B"
)

# Plot
ggplot(long_data, aes(x = N, y = Parameters, color = Case, shape = Case)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  labs(x = "Number of States", y = "Number of Parameters", color = "", shape = "") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
    legend.spacing.y = unit(1, "pt"),
    legend.text = element_text(margin = margin(t = -2))
  )
