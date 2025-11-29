# -----------------------------------------------
# Black–Scholes SSM: Viterbi + "like-their-analysis" plots (DateCont)
# -----------------------------------------------
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(grid)  # for unit()

# --- Inputs you already have ---
# spx             # data frame with columns: Return, DateCont
spx <- spx %>% slice_head(n = 23109)
 S <- na.omit(spx$Return)
x <- as.numeric(S)

# Align DateCont to the NA-omitted returns exactly like S
spx2 <- spx %>% slice_head(n = 23109)
idx_valid <- which(!is.na(spx2$Return))
DateCont  <- spx2$DateCont[idx_valid]
stopifnot(length(DateCont) == length(x))

# Parameters
theta_BS_hat <- c(phi   = 0.98416755,
                  sigZ  = 0.09326986,
                  mu    = 0.13905276,
                  sigma = 0.13140497)

dt <- 1/252
m  <- 200
bm <- 3

# Unpack
phi_state   <- as.numeric(theta_BS_hat["phi"])
sigma_state <- as.numeric(theta_BS_hat["sigZ"])
mu_hat      <- as.numeric(theta_BS_hat["mu"])
sigma_hat   <- as.numeric(theta_BS_hat["sigma"])

# Latent grid
b     <- seq(-bm, bm, length = m + 1)
h     <- b[2] - b[1]
bstar <- 0.5 * (b[-1] + b[-(m + 1)])
expZ  <- exp(bstar)

# Transition
Gamma <- matrix(0, m, m)
for (i in 1:m) {
  Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
}
Gamma <- sweep(Gamma, 1, pmax(rowSums(Gamma), 1e-300), "/")

# Stationary init
sd0   <- sigma_state / sqrt(1 - phi_state^2)
delta <- h * dnorm(bstar, mean = 0, sd = sd0)
delta <- delta / sum(delta)

# Emission per state
state_mean <- (mu_hat * expZ - 0.5 * (sigma_hat * expZ)^2) * dt
state_sd   <-  (sigma_hat * expZ) * sqrt(dt)

# ---------- Forward (for mixture weights) ----------
Tn   <- length(x)
phi  <- delta
wbar <- rep(0, m)

dens1 <- dnorm(x[1], mean = state_mean, sd = state_sd)
foo   <- phi * dens1
s1    <- sum(foo); stopifnot(is.finite(s1) && s1 > 0)
phi   <- foo / s1
wbar  <- wbar + phi

for (t in 2:Tn) {
  denst <- dnorm(x[t], mean = state_mean, sd = state_sd)
  foo   <- as.numeric(phi %*% Gamma) * denst
  st    <- sum(foo); if (!is.finite(st) || st <= 0) st <- 1e-300
  phi   <- foo / st
  wbar  <- wbar + phi
}
wbar <- wbar / Tn

# ---------- Viterbi ----------
logG   <- log(Gamma + 1e-300)
logpi0 <- log(delta + 1e-300)
logE   <- sapply(1:m, function(i) dnorm(x, mean = state_mean[i], sd = state_sd[i], log = TRUE))

V   <- matrix(-Inf, nrow = Tn, ncol = m)
psi <- matrix(0L,   nrow = Tn, ncol = m)
V[1, ] <- logpi0 + logE[1, ]

for (t in 2:Tn) {
  prev <- V[t - 1, ]
  for (j in 1:m) {
    vals        <- prev + logG[, j]
    i_star      <- which.max(vals)
    V[t, j]     <- vals[i_star] + logE[t, j]
    psi[t, j]   <- i_star
  }
}
path        <- integer(Tn)
path[Tn]    <- which.max(V[Tn, ])
for (t in (Tn - 1):1) path[t] <- psi[t + 1, path[t + 1]]

Z_hat      <- bstar[path]
expZ_hat   <- exp(Z_hat)
local_vol  <- sigma_hat * expZ_hat
local_mu   <- mu_hat   * expZ_hat
cond_mean  <- (local_mu - 0.5 * local_vol^2) * dt

# ---------- Data for plotting ----------
df <- tibble::tibble(
  DateCont   = as.numeric(DateCont),
  r          = x,
  Zhat       = Z_hat,
  expZhat    = expZ_hat,
  local_vol  = local_vol,
  cond_mean  = cond_mean
)

# Nice x breaks like in their plots (every ~5 years over range)
xmin <- floor(min(df$DateCont, na.rm = TRUE))
xmax <- ceiling(max(df$DateCont, na.rm = TRUE))
by   <- max(1, floor((xmax - xmin)/7))   # ~8 ticks like 1985,...,2020
x_breaks <- seq(xmin, xmax, by = by)




# ---- Plot 1: put colorbar on the RIGHT (outside) ----
pA <- ggplot(df, aes(x = DateCont, y = r, color = local_vol)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(x = "Time", y = "Log returns", color = expression(sigma*e^{Z[t]}),
       title = "Black–Scholes SSM w/ Log Returns") +
  scale_x_continuous(breaks = x_breaks) +
  scale_color_gradientn(colours = brewer.pal(n = 8, name = "Pastel2")) +
  guides(color = guide_colourbar(
    barwidth  = unit(0.25, "inches"),
    barheight = unit(3, "inches"),
    title.position = "top",
    ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE
  )) +
  theme(
    plot.title = element_text(size = 47, hjust = 0),
    axis.title = element_text(size = 37),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    legend.position = "right",          # <-- key change
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 16),
    plot.margin = margin(5.5, 20, 5.5, 5.5, "pt")  # extra right margin
  ) +
  coord_cartesian(clip = "off")  # ensures legend/margins display cleanly

print(pA)




# -----------------------------------------------
# Plot 3: Histogram + state-mixture density (time-averaged weights), same aesthetics
# -----------------------------------------------
xgrid <- seq(min(df$r), max(df$r), length.out = 1000)
mix_each <- sapply(1:m, function(i) wbar[i] * dnorm(xgrid, mean = state_mean[i], sd = state_sd[i]))
mix_tot  <- rowSums(mix_each)

# Thin component display for cleanliness
thin_idx  <- unique(round(seq(1, m, length.out = 30)))
mycolors  <- colorRampPalette(brewer.pal(n = 8, name = "Set3"))(length(thin_idx))

mix_df_total <- data.frame(x = xgrid, dens = mix_tot)
mix_df_each  <- do.call(rbind, lapply(seq_along(thin_idx), function(k){
  i <- thin_idx[k]
  data.frame(x = xgrid, dens = mix_each[, i], comp = paste0("S", i), ord = k)
}))

ggplot(df, aes(x = r, y = after_stat(density))) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6", alpha=0.7) +
  theme_bw() +
  xlab("Log returns") +
  ylab("Density") +
  ggtitle("Black–Scholes SSM w/ Log Returns") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, max(mix_tot) * 1.15)) +
  theme(plot.title = element_text(size=47, hjust=0),
        axis.title = element_text(size=37),
        axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = mix_df_each, aes(x = x, y = dens, group = comp),
            alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
  geom_line(data = mix_df_total, aes(x = x, y = dens),
            color = "black", linewidth = 1.3, linetype = "dotdash")



# Points colored by local drift mu * e^{Z_t}; legend on right
ggplot(df, aes(x = DateCont, y = r, color = local_mu)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(x = "Time", y = "Log returns", color = expression(mu*e^{Z[t]}),
       title = "Black–Scholes SSM: Returns colored by local drift") +
  scale_x_continuous(breaks = x_breaks) +
  # Drift can be +/- so use a diverging scale around 0
  scale_color_gradient2(midpoint = 0) +
  guides(color = guide_colourbar(
    barwidth = unit(0.25, "inches"),
    barheight = unit(3, "inches"),
    title.position = "top",
    ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE
  )) +
  theme(
    plot.title = element_text(size=47, hjust=0),
    axis.title = element_text(size=37),
    axis.text.x = element_text(size=23),
    axis.text.y = element_text(size=23),
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 16),
    plot.margin = margin(5.5, 20, 5.5, 5.5, "pt")
  )
