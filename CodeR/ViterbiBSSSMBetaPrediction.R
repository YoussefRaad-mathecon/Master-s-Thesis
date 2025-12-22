# BS–SSM_beta (AR(1) latent factor): OUT-OF-SAMPLE Viterbi + Decoded Plots + Histogram

# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(tidyverse)
library(ggplot2)
library(gridExtra)

# ---- Setup ----
set.seed(112233)

dt      <- 1/252
N_train <- 23109
m_grid  <- 100
bm_grid <- 4

# ---- Load fits ----
load("fitted_params_BSSSM_2019_beta.RData")  # expects: models_BSSSM_2019_beta

# ---- Load full returns ----
load("sp500_clean.Rdata")                    # expects: spx with Return and ideally Date

stopifnot(exists("spx"))
stopifnot("Return" %in% names(spx))

# ---- Out-of-sample sample ----
spx_full <- spx %>% filter(!is.na(Return))
if (nrow(spx_full) <= N_train) {
  stop("Not enough observations for out-of-sample: nrow(spx_full) <= N_train.")
}

spx_oos <- spx_full %>% slice((N_train + 1):nrow(spx_full))
x_ret   <- as.numeric(spx_oos$Return)

# ---- Time axis handling ----
spx_oos <- spx_oos %>% mutate(t_index = seq_len(n()))

use_date_axis <- ("Date" %in% names(spx_oos)) && inherits(spx_oos$Date, "Date")
x_col         <- if (use_date_axis) "Date" else "t_index"

scale_x_time <- function(df, x_col) {
  if (x_col == "Date") {
    scale_x_date(date_breaks = "2 years", date_labels = "%Y")
  } else {
    scale_x_continuous(breaks = pretty(df[[x_col]], n = 6))
  }
}

# ---- Reconstruct theta_BS_hat_beta ----
if (!exists("models_BSSSM_2019_beta")) {
  stop("models_BSSSM_2019_beta not found after loading fitted_params_BSSSM_2019_beta.RData")
}

if (!exists("theta_BS_hat_beta")) {
  fit_BS_beta <- NULL
  
  if (!is.null(models_BSSSM_2019_beta$fit_BS_beta)) {
    fit_BS_beta <- models_BSSSM_2019_beta$fit_BS_beta
  } else if (is.list(models_BSSSM_2019_beta) && length(models_BSSSM_2019_beta) >= 1) {
    fit_BS_beta <- models_BSSSM_2019_beta[[1]]
  }
  
  if (is.null(fit_BS_beta) || is.null(fit_BS_beta$estimate)) {
    stop("Could not recover nlm fit from models_BSSSM_2019_beta. Check the RData file structure.")
  }
  
  theta_hat_tr <- fit_BS_beta$estimate
  theta_BS_hat_beta <- c(
    phi        = plogis(theta_hat_tr[1]),
    sigZ       = exp(theta_hat_tr[2]),
    mu0        = theta_hat_tr[3],
    sigma      = exp(theta_hat_tr[4]),
    beta_mu    = theta_hat_tr[5],
    beta_sigma = theta_hat_tr[6]
  )
}

if (is.null(names(theta_BS_hat_beta))) {
  names(theta_BS_hat_beta) <- c("phi", "sigZ", "mu0", "sigma", "beta_mu", "beta_sigma")
}

# ---- Viterbi decoding ----
viterbi_bs_beta <- function(x, mu0, sigma, beta_mu, beta_sigma,
                            phi_state, sigma_state,
                            m = 100, bm = 4, dt = 1/252) {
  n <- length(x)
  
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  rs <- rowSums(Gamma)
  if (any(rs <= 0 | !is.finite(rs))) stop("Invalid Gamma rows in Viterbi (beta model).")
  Gamma <- Gamma / rs
  
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  sdel  <- sum(delta)
  if (!is.finite(sdel) || sdel <= 0) stop("Invalid stationary initial distribution in Viterbi.")
  delta <- delta / sdel
  
  mu_grid    <- mu0 + beta_mu * bstar
  sigma_grid <- sigma * exp(beta_sigma * bstar)
  
  mean_vec <- (mu_grid - 0.5 * sigma_grid^2) * dt
  sd_vec   <- sigma_grid * sqrt(dt)
  
  log_emiss <- function(x_t) dnorm(x_t, mean = mean_vec, sd = sd_vec, log = TRUE)
  
  logGamma <- log(Gamma)
  logdelta <- log(delta)
  
  nu  <- matrix(-Inf, n, m)
  psi <- matrix(0L, n, m)
  
  nu[1, ] <- logdelta + log_emiss(x[1])
  
  for (t in 2:n) {
    le <- log_emiss(x[t])
    for (j in 1:m) {
      v    <- nu[t - 1, ] + logGamma[, j]
      imax <- which.max(v)
      nu[t, j]  <- v[imax] + le[j]
      psi[t, j] <- imax
    }
  }
  
  iv <- integer(n)
  iv[n] <- which.max(nu[n, ])
  for (t in (n - 1):1) {
    iv[t] <- psi[t + 1, iv[t + 1]]
  }
  
  bstar[iv]
}

# ---- Out-of-sample decoding ----
C_hat_oos <- viterbi_bs_beta(
  x           = x_ret,
  mu0         = theta_BS_hat_beta["mu0"],
  sigma       = theta_BS_hat_beta["sigma"],
  beta_mu     = theta_BS_hat_beta["beta_mu"],
  beta_sigma  = theta_BS_hat_beta["beta_sigma"],
  phi_state   = theta_BS_hat_beta["phi"],
  sigma_state = theta_BS_hat_beta["sigZ"],
  m           = m_grid,
  bm          = bm_grid,
  dt          = dt
)

sigma_t_oos <- as.numeric(theta_BS_hat_beta["sigma"] * exp(theta_BS_hat_beta["beta_sigma"] * C_hat_oos))
mu_t_oos    <- as.numeric(theta_BS_hat_beta["mu0"]  + theta_BS_hat_beta["beta_mu"]    * C_hat_oos)

bs_beta_oos_df <- spx_oos %>%
  slice_head(n = length(C_hat_oos)) %>%
  mutate(
    C_hat   = C_hat_oos,
    sigma_t = sigma_t_oos,
    mu_t    = mu_t_oos
  )

# ---- Style ----
rainbow_palette <- grDevices::colorRampPalette(grDevices::rainbow(7))(256)

title_size      <- 20
axis_title_size <- 16
axis_text_size  <- 12

# ---- Plot: returns colored by sigma_t ----
sigma_range  <- range(bs_beta_oos_df$sigma_t, na.rm = TRUE)
sigma_breaks <- pretty(sigma_range, n = 5)

p_ret_sigma_oos <- ggplot(bs_beta_oos_df, aes(x = .data[[x_col]], y = Return, color = sigma_t)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(sigma[t]),
    limits  = sigma_range,
    breaks  = sigma_breaks
  ) +
  guides(
    colour = guide_colourbar(
      barheight = grid::unit(0.6, "npc"),
      barwidth  = grid::unit(0.02, "npc")
    )
  ) +
  theme_bw() +
  labs(
    title = expression(paste("BS-SSM" [beta], ": Volatility Decoded (Out-of-sample)")),
    x     = "Time",
    y     = "Daily log-Returns"
  ) +
  theme(
    plot.title      = element_text(size = title_size, hjust = 0),
    axis.title      = element_text(size = axis_title_size),
    axis.text       = element_text(size = axis_text_size),
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5, unit = "pt")
  ) +
  scale_x_time(bs_beta_oos_df, x_col)

p_ret_sigma_oos

# ---- Plot: returns colored by mu_t ----
mu_range  <- range(bs_beta_oos_df$mu_t, na.rm = TRUE)
mu_breaks <- pretty(mu_range, n = 5)

p_ret_mu_oos <- ggplot(bs_beta_oos_df, aes(x = .data[[x_col]], y = Return, color = mu_t)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(mu[t]),
    limits  = mu_range,
    breaks  = mu_breaks
  ) +
  guides(
    colour = guide_colourbar(
      barheight = grid::unit(0.6, "npc"),
      barwidth  = grid::unit(0.02, "npc")
    )
  ) +
  theme_bw() +
  labs(
    title = expression(paste("BS-SSM" [beta], ": Capital Drift Decoded (Out-of-sample)")),
    x     = "Time",
    y     = "Daily log-Returns"
  ) +
  theme(
    plot.title      = element_text(size = title_size, hjust = 0),
    axis.title      = element_text(size = axis_title_size),
    axis.text       = element_text(size = axis_text_size),
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5, unit = "pt")
  ) +
  scale_x_time(bs_beta_oos_df, x_col)

p_ret_mu_oos

# ---- Plot: latent factor C_t ----
Ct_range  <- range(bs_beta_oos_df$C_hat, na.rm = TRUE)
Ct_breaks <- pretty(Ct_range, n = 5)

p_Ct_oos <- ggplot(bs_beta_oos_df, aes(x = .data[[x_col]], y = C_hat, color = C_hat)) +
  geom_line(alpha = 0.9, linewidth = 0.4) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(C[t]),
    limits  = Ct_range,
    breaks  = Ct_breaks
  ) +
  guides(
    colour = guide_colourbar(
      barheight = grid::unit(0.6, "npc"),
      barwidth  = grid::unit(0.02, "npc")
    )
  ) +
  theme_bw() +
  labs(
    title = expression(paste("BS-SSM" [beta], ": Latent Factor Decoded (Out-of-sample)")),
    x     = "Time",
    y     = expression(C[t])
  ) +
  theme(
    plot.title      = element_text(size = title_size, hjust = 0),
    axis.title      = element_text(size = axis_title_size),
    axis.text       = element_text(size = axis_text_size),
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5, unit = "pt")
  ) +
  scale_x_time(bs_beta_oos_df, x_col) +
  scale_y_continuous(breaks = Ct_breaks)

p_Ct_oos


# ---- Save plots ----
ggsave("BSSSMbeta_OOS_VolatilityDecoded.png", plot = p_ret_sigma_oos, width = 10, height = 6, dpi = 300)
ggsave("BSSSMbeta_OOS_CapitalDriftDecoded.png", plot = p_ret_mu_oos,  width = 10, height = 6, dpi = 300)
ggsave("BSSSMbeta_OOS_LatentFactor.png",        plot = p_Ct_oos,      width = 10, height = 6, dpi = 300)
