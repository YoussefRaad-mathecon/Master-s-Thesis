# Viterbi + Plots: Black–Scholes BS-SSM_beta

# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

FIT_FILE  <- "fitted_params_BSSSM_2019_beta.RData"
DIV_FILE  <- "sp500_clean_with_divs.RData"
RET_FILE  <- "sp500_clean.Rdata"

set.seed(112233)

# ---- Packages ----
library(dplyr)
library(tibble)
library(ggplot2)
library(sandwich)

# ---- Load fitted SSM_beta object ----
load(FIT_FILE)  # expects: models_BSSSM_2019_beta (structure may vary)

# ---- Recover natural-scale parameters: (phi, sigZ, mu0, sigma, beta_mu, beta_sigma) ----
if (!exists("theta_BS_hat_beta")) {
  stopifnot(exists("models_BSSSM_2019_beta"))
  
  fit_BS_beta <- NULL
  
  if (!is.null(models_BSSSM_2019_beta$fit_BS_beta)) {
    fit_BS_beta <- models_BSSSM_2019_beta$fit_BS_beta
  } else if (!is.null(models_BSSSM_2019_beta$fit_object)) {
    fit_BS_beta <- models_BSSSM_2019_beta$fit_object
  } else if (is.list(models_BSSSM_2019_beta) && length(models_BSSSM_2019_beta) >= 1) {
    fit_BS_beta <- models_BSSSM_2019_beta[[1]]
  }
  
  if (is.null(fit_BS_beta) || is.null(fit_BS_beta$estimate)) {
    stop("Could not recover nlm fit from models_BSSSM_2019_beta. Check the RData structure.")
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

dt      <- 1/252
m_grid  <- 100
bm_grid <- 4

# ---- Load data (returns + dividends) ----
load(DIV_FILE)  # expects: spx_div with column 'DivLog'
load(RET_FILE)  # expects: spx with column 'Return' (and optionally Date, DateCont)

N_target <- 23109
N <- min(N_target, nrow(spx), nrow(spx_div))

spx     <- spx     %>% slice_head(n = N)
spx_div <- spx_div %>% slice_head(n = N)

valid_idx <- which(!is.na(spx$Return) & !is.na(spx_div$DivLog))

spx     <- spx[valid_idx, , drop = FALSE]
spx_div <- spx_div[valid_idx, , drop = FALSE]

x       <- as.numeric(spx$Return)
q_daily <- spx_div$DivLog

date_vec <- if ("Date" %in% names(spx)) spx$Date else seq_along(x)
year_vec <- if ("DateCont" %in% names(spx)) spx$DateCont else seq_along(x)

# ---- Constant annual dividend yield q_hat (HAC) ----
estimate_q_hac <- function(q_vec) {
  mod <- lm(q_vec ~ 1)
  bw  <- sandwich::bwNeweyWest(mod)
  vc  <- sandwich::NeweyWest(mod, lag = bw, prewhite = FALSE, adjust = TRUE)
  
  q_hat_daily <- coef(mod)[1]
  se_daily    <- sqrt(diag(vc))[1]
  
  list(
    q_hat_annual = 252 * q_hat_daily,
    se_annual    = 252 * se_daily
  )
}

q_est        <- estimate_q_hac(q_daily)
q_hat_annual <- q_est$q_hat_annual

# ---- Viterbi (grid HMM approximation of latent AR(1)) ----
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
  if (any(rs <= 0 | !is.finite(rs))) stop("Invalid Gamma rows in Viterbi.")
  Gamma <- Gamma / rs
  
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, 0, sd0)
  sdel  <- sum(delta)
  if (!is.finite(sdel) || sdel <= 0) stop("Invalid delta in Viterbi.")
  delta <- delta / sdel
  
  mu_grid    <- mu0 + beta_mu * bstar
  sigma_grid <- sigma * exp(beta_sigma * bstar)
  
  mean_vec <- (mu_grid - 0.5 * sigma_grid^2) * dt
  sd_vec   <- sigma_grid * sqrt(dt)
  
  log_emiss <- function(x_t) dnorm(x_t, mean = mean_vec, sd = sd_vec, log = TRUE)
  
  logGamma <- log(Gamma)
  logdelta <- log(delta)
  
  nu  <- matrix(-Inf, n, m)
  psi <- matrix(0L,   n, m)
  
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

# ---- Run Viterbi and build decoded paths ----
Ct <- viterbi_bs_beta(
  x           = x,
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

sigma_t <- as.numeric(theta_BS_hat_beta["sigma"] * exp(theta_BS_hat_beta["beta_sigma"] * Ct))
mu_t    <- as.numeric(theta_BS_hat_beta["mu0"]  + theta_BS_hat_beta["beta_mu"]    * Ct)
mu_tot_t <- mu_t + q_hat_annual

cond_mean_t <- (mu_t - 0.5 * sigma_t^2) * dt
cond_sd_t   <- sigma_t * sqrt(dt)

bs_beta_paths_df <- tibble(
  t           = seq_along(x),
  date        = date_vec,
  year        = year_vec,
  x           = x,
  C_hat       = Ct,
  sigma_t     = sigma_t,
  mu_t        = mu_t,
  mu_tot_t    = mu_tot_t,
  cond_mean_t = cond_mean_t,
  cond_sd_t   = cond_sd_t
)

# ---- Plots ----
rainbow_palette <- grDevices::colorRampPalette(grDevices::rainbow(7))(256)

x_breaks <- pretty(range(bs_beta_paths_df$year, na.rm = TRUE), n = 6)
y_breaks <- pretty(bs_beta_paths_df$x, n = 5)

title_front_sigma_beta <- expression(paste("BS-SSM" [beta], ": Volatility Decoded"))
title_front_mu_beta    <- expression(paste("BS-SSM" [beta], ": Capital Drift Decoded"))
title_front_Ct_beta    <- expression(paste("BS-SSM" [beta], ": Latent Factor Decoded"))

sigma_range  <- range(bs_beta_paths_df$sigma_t, na.rm = TRUE)
sigma_breaks <- pretty(sigma_range, n = 5)

p_front_sigma_beta <- ggplot(bs_beta_paths_df, aes(x = year, y = x, color = sigma_t)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(sigma[t]),
    limits  = sigma_range,
    breaks  = sigma_breaks
  ) +
  scale_x_continuous(breaks = x_breaks, labels = function(z) floor(z)) +
  scale_y_continuous(breaks = y_breaks) +
  guides(colour = guide_colourbar(
    barheight = grid::unit(0.6, "npc"),
    barwidth  = grid::unit(0.02, "npc")
  )) +
  theme_bw() +
  labs(title = title_front_sigma_beta, x = "Year", y = "Daily log-returns") +
  theme(
    plot.title      = element_text(size = 20, hjust = 0),
    panel.grid      = element_blank(),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  )

mu_range  <- range(bs_beta_paths_df$mu_t, na.rm = TRUE)
mu_breaks <- pretty(mu_range, n = 5)

p_front_mu_beta <- ggplot(bs_beta_paths_df, aes(x = year, y = x, color = mu_t)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(mu[t]),
    limits  = mu_range,
    breaks  = mu_breaks
  ) +
  scale_x_continuous(breaks = x_breaks, labels = function(z) floor(z)) +
  scale_y_continuous(breaks = y_breaks) +
  guides(colour = guide_colourbar(
    barheight = grid::unit(0.6, "npc"),
    barwidth  = grid::unit(0.02, "npc")
  )) +
  theme_bw() +
  labs(title = title_front_mu_beta, x = "Year", y = "Daily log-returns") +
  theme(
    plot.title      = element_text(size = 20, hjust = 0),
    panel.grid      = element_blank(),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  )

Ct_range  <- range(bs_beta_paths_df$C_hat, na.rm = TRUE)
Ct_breaks <- pretty(Ct_range, n = 5)

p_Ct_full_beta <- ggplot(bs_beta_paths_df, aes(x = year, y = C_hat, color = C_hat)) +
  geom_line(alpha = 0.9, linewidth = 0.4) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = expression(C[t]),
    limits  = Ct_range,
    breaks  = Ct_breaks
  ) +
  scale_x_continuous(breaks = x_breaks, labels = function(z) floor(z)) +
  scale_y_continuous(breaks = Ct_breaks) +
  guides(colour = guide_colourbar(
    barheight = grid::unit(0.6, "npc"),
    barwidth  = grid::unit(0.02, "npc")
  )) +
  theme_bw() +
  labs(title = title_front_Ct_beta, x = "Year", y = expression(C[t])) +
  theme(
    plot.title      = element_text(size = 20, hjust = 0),
    panel.grid      = element_blank(),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  )

# ---- Summary printouts ----
cat("\n================ BS-SSM_beta: Summary statistics for mu_t and sigma_t ================\n")

cat("\n--- mu_t (capital-gains drift) ---\n")
print(summary(bs_beta_paths_df$mu_t))
idx_mu_min <- which.min(bs_beta_paths_df$mu_t)
idx_mu_max <- which.max(bs_beta_paths_df$mu_t)
cat(sprintf("Min(mu_t) = %.4f at t = %d, year = %s, date = %s\n",
            bs_beta_paths_df$mu_t[idx_mu_min], idx_mu_min,
            as.character(bs_beta_paths_df$year[idx_mu_min]),
            as.character(bs_beta_paths_df$date[idx_mu_min])))
cat(sprintf("Max(mu_t) = %.4f at t = %d, year = %s, date = %s\n",
            bs_beta_paths_df$mu_t[idx_mu_max], idx_mu_max,
            as.character(bs_beta_paths_df$year[idx_mu_max]),
            as.character(bs_beta_paths_df$date[idx_mu_max])))

cat("\n--- sigma_t (volatility) ---\n")
print(summary(bs_beta_paths_df$sigma_t))
idx_sig_min <- which.min(bs_beta_paths_df$sigma_t)
idx_sig_max <- which.max(bs_beta_paths_df$sigma_t)
cat(sprintf("Min(sigma_t) = %.4f at t = %d, year = %s, date = %s\n",
            bs_beta_paths_df$sigma_t[idx_sig_min], idx_sig_min,
            as.character(bs_beta_paths_df$year[idx_sig_min]),
            as.character(bs_beta_paths_df$date[idx_sig_min])))
cat(sprintf("Max(sigma_t) = %.4f at t = %d, year = %s, date = %s\n",
            bs_beta_paths_df$sigma_t[idx_sig_max], idx_sig_max,
            as.character(bs_beta_paths_df$year[idx_sig_max]),
            as.character(bs_beta_paths_df$date[idx_sig_max])))

mu_cap_annual_beta <- as.numeric(theta_BS_hat_beta["mu0"])
mu_tot_annual_beta <- mu_cap_annual_beta + q_hat_annual
cat("BS-SSM_beta: mu_tot (annual, at C_t = 0) =", sprintf("%.4f", mu_tot_annual_beta), "\n")

# ---- Minimalist front-page C_t ----
p_front_Ct_beta <- ggplot(bs_beta_paths_df, aes(x = year, y = C_hat, color = C_hat)) +
  geom_line(alpha = 0.9, linewidth = 0.4) +
  scale_color_gradientn(
    colours = rainbow_palette,
    name    = NULL,
    limits  = Ct_range,
    breaks  = NULL
  ) +
  scale_x_continuous(breaks = x_breaks, labels = function(z) floor(z)) +
  scale_y_continuous(breaks = Ct_breaks) +
  guides(colour = guide_colourbar(
    barheight = grid::unit(0.9, "npc"),
    barwidth  = grid::unit(0.02, "npc")
  )) +
  theme_bw() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    panel.grid      = element_blank(),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    legend.position = "right",
    plot.margin     = margin(5, 5, 5, 5)
  )

# ---- Optional saves ----
# ggsave("BSSSM_beta_VolatilityDecoded.png",     plot = p_front_sigma_beta, width = 10, height = 6, dpi = 300)
# ggsave("BSSSM_beta_CapitalDriftDecoded.png",   plot = p_front_mu_beta,    width = 10, height = 6, dpi = 300)
# ggsave("BSSSM_beta_LatentFactorDecoded.png",   plot = p_Ct_full_beta,     width = 10, height = 6, dpi = 300)
# ggsave("BSSSM_beta_Ct_Frontpage.png",          plot = p_front_Ct_beta,    width = 10, height = 6, dpi = 300)
