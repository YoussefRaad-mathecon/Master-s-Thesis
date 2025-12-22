# BS–HMM (4 states, μ & σ): OUT-OF-SAMPLE Viterbi + State Path + Histogram

library(tidyverse)
library(matrixStats)
library(ggplot2)
library(gridExtra)

set.seed(112233)

###--------------------------------------------------------------------------------------------------
### 0) Paths & OUT-OF-SAMPLE data
###--------------------------------------------------------------------------------------------------
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("fitted_params_BSHMM_2019.RData")  # list: models_BSHMM_2019
load("sp500_clean.Rdata")               # expects `spx$Return` and `spx$Date`

# Full sample, NA-free in Return
spx_full <- spx %>% filter(!is.na(Return))

# Out-of-sample (everything after the 23109 obs used for estimation)
spx <- spx_full %>% slice((23109 + 1):n())

x_ret <- as.numeric(spx$Return)
Delta <- 1/252

# For histogram section we also use a clean, NA-free vector (OOS only)
x  <- as.numeric(na.omit(spx$Return))
dt <- 1/252

###--------------------------------------------------------------------------------------------------
### 1) Style (shared)
###--------------------------------------------------------------------------------------------------
# 5-color palette (first 4 used for 4-state BS-HMM)
palette_states <- c("#DEAB7E", "#78C2AF", "#77ACC9", "#D98E96", "#B39AD8")

################################################################################################
# 2) Viterbi decoding for 4-state BS–HMM (μ & σ) on OOS
################################################################################################

viterbi <- function(x, mu, sigma, Gamma, delta, N){
  dt <- 1/252
  n  <- length(x)
  allprobs <- matrix(1, nrow = n, ncol = N)
  ind <- which(!is.na(x))
  
  for (j in 1:N){
    allprobs[ind, j] <- dnorm(
      x[ind],
      mean = (mu[j] - 0.5 * sigma[j]^2) * dt,
      sd   = sqrt(dt * sigma[j]^2)
    )
  }
  
  xi  <- matrix(0, nrow = n, ncol = N)
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  
  for (t in 2:n){
    foo      <- apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ]  <- foo / sum(foo)
  }
  
  iv     <- numeric(n)
  iv[n]  <- which.max(xi[n, ])
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ])
  }
  iv
}

### Helper: unpack a BS–HMM parameter vector
unpack_bshmm <- function(fit, N, family = c("mu","sigma","mu_sigma")) {
  family <- match.arg(family)
  th     <- fit$estimate
  
  if (family == "mu") {
    mu        <- th[1:N]
    log_sigma <- th[N + 1]
    eta       <- th[(N + 2):length(th)]
    sigma     <- rep(exp(log_sigma), N)
  } else if (family == "sigma") {
    mu        <- rep(th[1], N)
    log_sigma <- th[2:(N + 1)]
    eta       <- th[(N + 2):length(th)]
    sigma     <- exp(log_sigma)
  } else { # "mu_sigma"
    mu        <- th[1:N]
    log_sigma <- th[(N + 1):(2 * N)]
    eta       <- th[(2 * N + 1):length(th)]
    sigma     <- exp(log_sigma)
  }
  
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- sweep(Gamma, 1, rowSums(Gamma), "/")
  
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  delta <- as.numeric(delta / sum(delta))
  
  list(mu = mu, sigma = sigma, Gamma = Gamma, delta = delta)
}

# Decode parameters for the 4-state μ & σ model
N             <- 4
mod4_mu_sigma <- models_BSHMM_2019$mod4_mu_sigma
pars4         <- unpack_bshmm(mod4_mu_sigma, N, family = "mu_sigma")

mu4_hat    <- pars4$mu
sigma4_hat <- pars4$sigma
Gamma4_hat <- pars4$Gamma
delta4_hat <- pars4$delta

# Viterbi decoding on OOS x_ret
viterbi4 <- viterbi(
  x     = x_ret,
  mu    = mu4_hat,
  sigma = sigma4_hat,
  Gamma = Gamma4_hat,
  delta = delta4_hat,
  N     = N
)

# Data frame with states (OOS)
spx_states <- spx %>%
  slice_head(n = length(viterbi4)) %>%
  mutate(Viterbi4 = viterbi4)

################################################################################################
# 2b) OUT-OF-SAMPLE legend numbers (decoded shares)
################################################################################################

v_oos <- as.integer(viterbi4)

state_pct_oos <- 100 * tabulate(v_oos, nbins = N) / length(v_oos)

# Round to 3 significant digits and force sum to 100
state_pct_oos_round <- signif(state_pct_oos, 3)
if (N > 1) {
  state_pct_oos_round[N] <- 100 - sum(state_pct_oos_round[1:(N - 1)])
}

state_labs_oos <- paste0(
  "S", seq_len(N), ": ",
  formatC(
    state_pct_oos_round,
    format        = "fg",
    digits        = 3,
    drop0trailing = TRUE,
    decimal.mark  = "."
  ),
  "%"
)

################################################################################################
# 3) Plots: BS–HMM State Decoding + State Path (OOS)
################################################################################################

## 3a) State decoding: returns colored by state
x_range_dates <- range(spx_states$Date, na.rm = TRUE)

# Put legend around 80% of the time axis (inside panel, away from axes)
x_legend <- x_range_dates[1] + 0.80 * diff(x_range_dates)

# y positioning: relative to data range so spacing looks ok
y_range      <- range(spx_states$Return, na.rm = TRUE)
y_span       <- diff(y_range)
y_legend_top <- y_range[2] - 0.05 * y_span
y_step       <- 0.06 * y_span

p_ret <- spx_states %>%
  ggplot(aes(x = Date, y = Return)) +
  geom_point(aes(color = factor(Viterbi4)), alpha = 0.5, size = 2) +
  theme_bw() +
  xlab("Time") +
  ylab("Daily log-Returns") +
  ggtitle("BS-HMM State Decoding (Out-of-sample)") +
  theme(
    plot.title   = element_text(size = 20, hjust = 0),
    axis.title   = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    plot.margin  = margin(5, 45, 5, 5, unit = "pt")
  ) +
  scale_color_manual(values = palette_states, guide = "none") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  annotate("text", x = x_legend, y = y_legend_top,             label = state_labs_oos[1], color = "#DEAB7E", size = 5, hjust = 0) +
  annotate("text", x = x_legend, y = y_legend_top - y_step,     label = state_labs_oos[2], color = "#78C2AF", size = 5, hjust = 0) +
  annotate("text", x = x_legend, y = y_legend_top - 2 * y_step, label = state_labs_oos[3], color = "#77ACC9", size = 5, hjust = 0) +
  annotate("text", x = x_legend, y = y_legend_top - 3 * y_step, label = state_labs_oos[4], color = "#D98E96", size = 5, hjust = 0)

p_ret

## 3b) State path: discrete states over time
x_legend_state <- x_range_dates[1] + 0.80 * diff(x_range_dates)
y_legend_top   <- 4.6
y_step         <- 0.18

p_state_path <- spx_states %>%
  ggplot(aes(x = Date, y = Viterbi4)) +
  geom_point(aes(color = factor(Viterbi4)), alpha = 0.5, size = 2) +
  theme_bw() +
  xlab("Time") +
  ylab("State") +
  ggtitle("BS-HMM State Path (Out-of-sample)") +
  theme(
    plot.title   = element_text(size = 20, hjust = 0),
    axis.title   = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    plot.margin  = margin(5, 45, 5, 5, unit = "pt")
  ) +
  scale_color_manual(values = palette_states, guide = "none") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(breaks = 1:4, limits = c(1, 4.8)) +
  annotate("text", x = x_legend_state, y = y_legend_top,             label = state_labs_oos[1], color = palette_states[1], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top - y_step,     label = state_labs_oos[2], color = palette_states[2], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top - 2 * y_step, label = state_labs_oos[3], color = palette_states[3], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top - 3 * y_step, label = state_labs_oos[4], color = palette_states[4], size = 5, hjust = 0)

p_state_path

################################################################################################
# 4) Histogram + State-dependent BS–HMM Return Densities + Mixture  (OOS)
################################################################################################

title_size      <- 18
axis_title_size <- 16
axis_text_size  <- 12
annot_size      <- 5
line_width      <- 0.9
hist_alpha      <- 0.55

# Global x-range & grid for densities (OOS)
xrng_all <- quantile(x, c(0.001, 0.999), na.rm = TRUE)
x_grid   <- seq(xrng_all[1], xrng_all[2], length.out = 1000)

title_for <- function(N, fam) {
  if (fam == "mu")        return(paste0(N, "-State BS-HMM (", "\u03BC", ") (Out-of-Sample)"))
  if (fam == "sigma")     return(paste0(N, "-State BS-HMM (", "\u03C3", ") (Out-of-Sample)"))
  if (fam == "mu_sigma")  return(paste0(N, "-State BS-HMM (", "\u03BC", " & ", "\u03C3", ") (Out-of-Sample)"))
  stop("Unknown family")
}

# ---- ONLY CHANGE HERE: add `legend_labs` and use it for labels_df ----
plot_state_densities <- function(x, fit, N, family, legend_labs) {
  pars  <- unpack_bshmm(fit, N, family)
  means <- (pars$mu - 0.5 * pars$sigma^2) * dt
  sds   <- pars$sigma * sqrt(dt)
  
  xrng   <- xrng_all
  xrange <- diff(xrng)
  
  x_label_text <- xrng[1] + 0.86 * xrange
  x_mask_min   <- x_label_text - 0.02 * xrange
  x_mask_max   <- xrng[2]
  
  xg_data <- x_grid
  
  dens_df <- purrr::map_dfr(seq_len(N), function(i) {
    tibble(
      x       = xg_data,
      density = dnorm(xg_data, mean = means[i], sd = sds[i]),
      i       = i
    )
  })
  
  mix_df <- tibble(
    x       = xg_data,
    density = rowSums(sapply(seq_len(N), function(i) {
      pars$delta[i] * dnorm(xg_data, mean = means[i], sd = sds[i])
    }))
  )
  
  h  <- hist(x, breaks = 60, plot = FALSE, include.lowest = TRUE, right = TRUE)
  bw <- diff(h$breaks)
  hist_densities    <- h$counts / sum(h$counts) / bw
  hist_density_max  <- max(hist_densities, na.rm = TRUE)
  
  line_density_max  <- max(dens_df$density, na.rm = TRUE)
  mix_density_max   <- max(mix_df$density,  na.rm = TRUE)
  y_max_base        <- max(line_density_max, mix_density_max, hist_density_max)
  
  if (N <= 2) {
    y_max_factor <- 1.15
  } else if (N == 3) {
    y_max_factor <- 1.20
  } else if (N == 4) {
    y_max_factor <- 1.25
  } else {
    y_max_factor <- 1.30
  }
  y_max_raw <- y_max_base * y_max_factor
  y_max     <- min(y_max_raw, 100)
  
  y_mask_min <- min(hist_density_max * 1.02, y_max)
  
  cols <- palette_states[seq_len(N)]
  
  top_frac <- 0.88
  gap_frac <- 0.048
  
  y_top_labels <- y_max * top_frac
  y_vals       <- y_top_labels - (seq_len(N) - 1L) * (gap_frac * y_max)
  
  # ---- FIXED: use the same legend labels as the other plots ----
  labels_df <- tibble(
    x   = x_label_text,
    y   = y_vals,
    lab = legend_labs,
    i   = seq_len(N)
  )
  
  header_y <- y_top_labels + 0.045 * y_max
  
  ggplot() +
    geom_histogram(
      data = data.frame(x = x),
      aes(x = x, y = after_stat(density)),
      bins  = 60,
      color = "#4D5D53",
      fill  = "#F0EAD6",
      alpha = hist_alpha
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = density, color = factor(i)),
      linewidth = line_width
    ) +
    geom_line(
      data = mix_df,
      aes(x = x, y = density),
      color    = "black",
      linewidth = line_width + 0.2,
      linetype = "dotted"
    ) +
    geom_rect(
      aes(xmin = x_mask_min, xmax = x_mask_max, ymin = y_mask_min, ymax = y_max),
      inherit.aes = FALSE,
      fill  = NA,
      color = NA
    ) +
    geom_text(
      data = labels_df,
      aes(x = x, y = y, label = lab, color = factor(i)),
      size  = annot_size,
      hjust = 0,
      vjust = 0.5
    ) +
    annotate(
      "text",
      x = x_label_text,
      y = header_y,
      label = "Fitted model",
      hjust  = 0,
      vjust  = 0.5,
      size   = annot_size,
      color  = "black"
    ) +
    scale_color_manual(values = cols, guide = "none") +
    theme_bw() +
    labs(
      title = title_for(N, family),
      x     = "Daily log-returns",
      y     = "Density"
    ) +
    scale_x_continuous(
      breaks = pretty(xrng_all, n = 5),
      limits = xrng_all,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 25),
      limits = c(0, 100),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title       = element_text(size = title_size, hjust = 0),
      axis.title       = element_text(size = axis_title_size),
      axis.text        = element_text(size = axis_text_size),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(5, 45, 5, 5, unit = "pt")
    )
}

# Build 4-state (μ & σ) histogram panel (OOS)  -- pass the OOS labels here
fit_4_mu_sigma    <- mod4_mu_sigma
p_hist_4_mu_sigma <- plot_state_densities(
  x           = x,
  fit         = fit_4_mu_sigma,
  N           = 4,
  family      = "mu_sigma",
  legend_labs = state_labs_oos
)

p_hist_4_mu_sigma

################################################################################################
# Save to disk
################################################################################################
ggsave("BSHMM_OOS_StateDecoding.png", plot = p_ret,             width = 10, height = 6, dpi = 300)
ggsave("BSHMM_OOS_StatePath.png",     plot = p_state_path,      width = 10, height = 6, dpi = 300)
ggsave("BSHMM_OOS_Histogram.png",     plot = p_hist_4_mu_sigma, width = 10, height = 6, dpi = 300)
