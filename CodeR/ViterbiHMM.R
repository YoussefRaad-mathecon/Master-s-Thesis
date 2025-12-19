# Viterbi + Plots: Black–Scholes BS-HMM


# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

set.seed(112233)

# ---- Load data ----
FIT_FILE  <- if (file.exists("fitted_params_BSHMM_2019.RData")) "fitted_params_BSHMM_2019.RData" else "fitted_params_BSHMM_2019.Rdata"
DATA_FILE <- if (file.exists("sp500_clean.Rdata")) "sp500_clean.Rdata" else "sp500_clean.RData"

stopifnot(file.exists(FIT_FILE), file.exists(DATA_FILE))

load(FIT_FILE)   # expects: models_BSHMM_2019
load(DATA_FILE)  # expects: spx with column Return (and ideally Date)

stopifnot(exists("models_BSHMM_2019"), exists("spx"))
stopifnot("Return" %in% names(spx))

N_target <- 23109
spx <- spx %>% slice_head(n = min(N_target, nrow(spx)))

dt <- 1/252

# Return vector used for Viterbi (keeps NAs to match Date alignment)
x_ret <- as.numeric(spx$Return)

# NA-free vector used for histogram density fitting
x <- as.numeric(na.omit(spx$Return))

# Ensure a Date axis exists
if (!("Date" %in% names(spx))) {
  spx$Date <- as.Date("1900-01-01") + seq_len(nrow(spx)) - 1L
}

# ---- Style ----
palette_states <- c("#DEAB7E", "#78C2AF", "#77ACC9", "#D98E96", "#B39AD8")

# ---- Viterbi decoding ----
viterbi <- function(x, mu, sigma, Gamma, delta, N, dt = 1/252) {
  n <- length(x)
  
  allprobs <- matrix(1, nrow = n, ncol = N)
  ind <- which(!is.na(x))
  
  for (j in 1:N) {
    allprobs[ind, j] <- dnorm(
      x[ind],
      mean = (mu[j] - 0.5 * sigma[j]^2) * dt,
      sd   = sigma[j] * sqrt(dt)
    )
  }
  
  xi  <- matrix(0, nrow = n, ncol = N)
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  
  for (t in 2:n) {
    foo     <- apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo / sum(foo)
  }
  
  iv    <- integer(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ])
  }
  iv
}

# ---- Unpack parameters ----
unpack_bshmm <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
  family <- match.arg(family)
  th <- fit$estimate
  
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
    
  } else {
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

# ---- Legend labels ----
state_labels_from_delta <- function(delta, digits = 3) {
  N <- length(delta)
  pct <- 100 * delta
  pct_round <- signif(pct, digits)
  if (N > 1) pct_round[N] <- 100 - sum(pct_round[1:(N - 1)])
  paste0(
    "S", seq_len(N), ": ",
    formatC(pct_round, format = "fg", digits = digits, drop0trailing = TRUE, decimal.mark = "."),
    "%"
  )
}

# ---- Decode 4-state (mu & sigma) model ----
N <- 4
stopifnot("mod4_mu_sigma" %in% names(models_BSHMM_2019))

fit_4_mu_sigma <- models_BSHMM_2019$mod4_mu_sigma
pars4 <- unpack_bshmm(fit_4_mu_sigma, N = N, family = "mu_sigma")

mu4_hat    <- pars4$mu
sigma4_hat <- pars4$sigma
Gamma4_hat <- pars4$Gamma
delta4_hat <- pars4$delta

viterbi4 <- viterbi(
  x     = x_ret,
  mu    = mu4_hat,
  sigma = sigma4_hat,
  Gamma = Gamma4_hat,
  delta = delta4_hat,
  N     = N,
  dt    = dt
)

spx_states <- spx %>%
  slice_head(n = length(viterbi4)) %>%
  mutate(Viterbi4 = viterbi4)

state_labs <- state_labels_from_delta(delta4_hat)

# ---- Plot: state decoding on returns ----
x_legend_ret <- max(spx_states$Date) - 7 * 365

y_rng <- range(spx_states$Return, na.rm = TRUE)
y_top <- y_rng[2]
y_span <- diff(y_rng)
y_legend_top <- y_top - 0.05 * y_span
y_step <- 0.06 * y_span

p_ret <- ggplot(spx_states, aes(x = Date, y = Return)) +
  geom_point(aes(color = factor(Viterbi4)), alpha = 0.5, size = 2) +
  theme_bw() +
  labs(
    title = "BS-HMM State Decoding",
    x = "Time",
    y = "Daily log-returns"
  ) +
  theme(
    plot.title  = element_text(size = 20, hjust = 0),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_manual(values = palette_states[1:N], guide = "none") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  annotate("text", x = x_legend_ret, y = y_legend_top - 0 * y_step, label = state_labs[1], color = palette_states[1], size = 5, hjust = 0) +
  annotate("text", x = x_legend_ret, y = y_legend_top - 1 * y_step, label = state_labs[2], color = palette_states[2], size = 5, hjust = 0) +
  annotate("text", x = x_legend_ret, y = y_legend_top - 2 * y_step, label = state_labs[3], color = palette_states[3], size = 5, hjust = 0) +
  annotate("text", x = x_legend_ret, y = y_legend_top - 3 * y_step, label = state_labs[4], color = palette_states[4], size = 5, hjust = 0)

p_ret

# ---- Plot: Viterbi state path ----
x_legend_state <- max(spx_states$Date) - 7 * 365
y_legend_top_state <- N + 0.6
y_step_state <- 0.18

p_state_path <- ggplot(spx_states, aes(x = Date, y = Viterbi4)) +
  geom_point(aes(color = factor(Viterbi4)), alpha = 0.5, size = 2) +
  theme_bw() +
  labs(
    title = "BS-HMM State Path",
    x = "Time",
    y = "State"
  ) +
  theme(
    plot.title  = element_text(size = 20, hjust = 0),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_manual(values = palette_states[1:N], guide = "none") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  scale_y_continuous(breaks = 1:N, limits = c(1, N + 0.8)) +
  annotate("text", x = x_legend_state, y = y_legend_top_state - 0 * y_step_state, label = state_labs[1], color = palette_states[1], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top_state - 1 * y_step_state, label = state_labs[2], color = palette_states[2], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top_state - 2 * y_step_state, label = state_labs[3], color = palette_states[3], size = 5, hjust = 0) +
  annotate("text", x = x_legend_state, y = y_legend_top_state - 3 * y_step_state, label = state_labs[4], color = palette_states[4], size = 5, hjust = 0)

p_state_path

# ---- Histogram + state densities + mixture ----
title_size      <- 18
axis_title_size <- 16
axis_text_size  <- 12
annot_size      <- 5
line_width      <- 0.9
hist_alpha      <- 0.55

xrng_all <- as.numeric(quantile(x, c(0.001, 0.999), na.rm = TRUE))
x_grid   <- seq(xrng_all[1], xrng_all[2], length.out = 1000)

title_for <- function(N, fam) {
  if (fam == "mu")       return(paste0(N, "-State BS-HMM (\u03BC)"))
  if (fam == "sigma")    return(paste0(N, "-State BS-HMM (\u03C3)"))
  if (fam == "mu_sigma") return(paste0(N, "-State BS-HMM (\u03BC & \u03C3)"))
  stop("Unknown family")
}

plot_state_densities <- function(x, fit, N, family, dt = 1/252) {
  pars  <- unpack_bshmm(fit, N, family)
  means <- (pars$mu - 0.5 * pars$sigma^2) * dt
  sds   <- pars$sigma * sqrt(dt)
  
  xrng   <- xrng_all
  xrange <- diff(xrng)
  
  x_label_text <- xrng[1] + 0.88 * xrange
  x_mask_min   <- x_label_text - 0.02 * xrange
  x_mask_max   <- xrng[2]
  
  dens_df <- purrr::map_dfr(seq_len(N), function(i) {
    tibble(
      x       = x_grid,
      density = dnorm(x_grid, mean = means[i], sd = sds[i]),
      i       = i
    )
  })
  
  mix_df <- tibble(
    x       = x_grid,
    density = rowSums(sapply(seq_len(N), function(i) {
      pars$delta[i] * dnorm(x_grid, mean = means[i], sd = sds[i])
    }))
  )
  
  h  <- hist(x, breaks = 60, plot = FALSE, include.lowest = TRUE, right = TRUE)
  bw <- diff(h$breaks)
  hist_densities   <- h$counts / sum(h$counts) / bw
  hist_density_max <- max(hist_densities, na.rm = TRUE)
  
  line_density_max <- max(dens_df$density, na.rm = TRUE)
  mix_density_max  <- max(mix_df$density,  na.rm = TRUE)
  y_max_base       <- max(line_density_max, mix_density_max, hist_density_max)
  
  y_max_factor <- if (N <= 2) 1.15 else if (N == 3) 1.20 else if (N == 4) 1.25 else 1.30
  y_max_raw <- y_max_base * y_max_factor
  y_max     <- min(y_max_raw, 100)
  
  y_mask_min <- min(hist_density_max * 1.02, y_max)
  
  cols <- palette_states[seq_len(N)]
  
  top_frac <- 0.88
  gap_frac <- 0.048
  
  y_top_labels <- y_max * top_frac
  y_vals       <- y_top_labels - (seq_len(N) - 1L) * (gap_frac * y_max)
  
  labels_df <- tibble(
    x   = x_label_text,
    y   = y_vals,
    lab = paste0("S", seq_len(N), ": ", sprintf("%.1f%%", 100 * pars$delta)),
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
      fill  = "white",
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
      hjust = 0,
      vjust = 0.5,
      size  = annot_size,
      color = "black"
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

p_hist_4_mu_sigma <- plot_state_densities(
  x      = x,
  fit    = fit_4_mu_sigma,
  N      = 4,
  family = "mu_sigma",
  dt     = dt
)

p_hist_4_mu_sigma
