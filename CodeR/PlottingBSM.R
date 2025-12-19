# BSM (No Regimes): Histogram + fitted Normal density

# ---- Packages ----
library(tidyverse)
library(ggplot2)

# ---- Data and Path ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("sp500_clean.Rdata")

x  <- as.numeric(na.omit(spx$Return))
dt <- 1 / 252

# ---- Black–Scholes parameters (annual) ----
mu_bs    <- 0.074541
sigma_bs <- 0.188326

# Daily mean and sd under BS
mean_bs <- (mu_bs - 0.5 * sigma_bs^2) * dt
sd_bs   <- sigma_bs * sqrt(dt)

# ---- Density grid + range ----
xrng_all <- quantile(x, c(0.001, 0.999), na.rm = TRUE)
x_grid   <- seq(xrng_all[1], xrng_all[2], length.out = 1000)

dens_bs <- tibble(
  x       = x_grid,
  density = dnorm(x_grid, mean = mean_bs, sd = sd_bs)
)

# ---- Style ----
hist_alpha      <- 0.55
line_width      <- 0.9
title_size      <- 18
axis_title_size <- 13
axis_text_size  <- 10

# ---- Plot: histogram + BS density ----
p_bs <- ggplot() +
  geom_histogram(
    data = data.frame(x = x),
    aes(x = x, y = after_stat(density)),
    bins  = 60,
    color = "#4D5D53",
    fill  = "#F0EAD6",
    alpha = hist_alpha
  ) +
  geom_line(
    data = dens_bs,
    aes(x = x, y = density),
    linewidth = line_width,
    color     = "black"
  ) +
  theme_bw() +
  labs(
    title = "BSM (No Regimes)",
    x     = "Log-returns",
    y     = "Density"
  ) +
  scale_x_continuous(
    breaks = pretty(xrng_all, n = 5),
    limits = xrng_all,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.02)
  ) +
  theme(
    plot.title       = element_text(size = title_size, hjust = 0),
    axis.title       = element_text(size = axis_title_size),
    axis.text        = element_text(size = axis_text_size),
    panel.grid.minor = element_blank()
  )

print(p_bs)
