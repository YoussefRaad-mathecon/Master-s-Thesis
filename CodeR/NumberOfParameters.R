# BS–HMM: Parameter-count plot by number of states (N)

# ---- Packages ----
library(ggplot2)
library(tidyr)
library(dplyr)

# ---- Setup ----
N <- 1:6

# ---- Parameter counts under your BS–HMM parameterization ----
# 2 state-dependent (μ_i and σ_i):  N(N−1) + 2N  = N^2 + N
# 1 state-dependent (exactly one):  N(N−1) + (N+1) = N^2 + 1
# 0 state-dependent (neither):      N(N−1) + 2   = N^2 − N + 2
data <- data.frame(
  N = N,
  Parameters_2 = N^2 + N,          # both μ and σ state-dependent
  Parameters_1 = N^2 + 1,          # exactly one of μ or σ state-dependent
  Parameters_0 = N^2 - N + 2       # neither state-dependent
)

# ---- Long format for ggplot ----
long_data <- pivot_longer(
  data,
  cols = starts_with("Parameters"),
  names_to = "Case",
  values_to = "Parameters"
)

# ---- Legend labels ----
long_data$Case <- recode(
  long_data$Case,
  Parameters_2 = "2 state-dependent",
  Parameters_1 = "1 state-dependent",
  Parameters_0 = "0 state-dependent"
)

# ---- Aesthetics ----
shape_map <- c(
  "2 state-dependent" = 15,  # square
  "1 state-dependent" = 17,  # triangle
  "0 state-dependent" = 16   # circle
)

color_map <- c(
  "2 state-dependent" = "#CC5500",
  "1 state-dependent" = "#404080",
  "0 state-dependent" = "#87A96B"
)

# ---- Plot ----
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
