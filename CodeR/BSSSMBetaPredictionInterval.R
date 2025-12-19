# BS-SSM_beta (Black–Scholes SSM with AR(1) latent factor): Forecast table + 1-step-ahead MSE/RMSE/MAE

# ---- Packages ----
library(tidyverse)
library(lubridate)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data, model, and train/test split ----
load("sp500_clean.Rdata")                     # expects spx$Date, spx$Return
load("fitted_params_BSSSM_2019_beta.RData")   # expects list models_BSSSM_2019_beta

if ("Date" %in% names(spx)) {
  spx <- spx %>% arrange(Date)
}

x_all <- as.numeric(spx$Return)
dates <- if ("Date" %in% names(spx)) spx$Date else seq_along(x_all)
n_all <- length(x_all)

n_train <- 23109L
if (n_all <= n_train) stop("Not enough observations for a test set: n_total <= 23109.")

idx_train <- 1:n_train
idx_test  <- (n_train + 1):n_all

x_train <- x_all[idx_train]
x_test  <- x_all[idx_test]

date_train <- dates[idx_train]
date_test  <- dates[idx_test]

# ---- Recover natural-scale parameters and discretisation settings ----
stopifnot(exists("models_BSSSM_2019_beta"))
stopifnot(!is.null(models_BSSSM_2019_beta$theta_hat))

theta_hat_nat <- models_BSSSM_2019_beta$theta_hat

# Support both named and unnamed theta_hat
get_th <- function(th, nm, pos) {
  if (!is.null(names(th)) && nm %in% names(th)) return(as.numeric(th[nm]))
  as.numeric(th[pos])
}

phi_hat        <- get_th(theta_hat_nat, "phi",        1)
sigZ_hat       <- get_th(theta_hat_nat, "sigZ",       2)
mu0_hat        <- get_th(theta_hat_nat, "mu0",        3)
sigma_hat      <- get_th(theta_hat_nat, "sigma",      4)
beta_mu_hat    <- get_th(theta_hat_nat, "beta_mu",    5)
beta_sigma_hat <- get_th(theta_hat_nat, "beta_sigma", 6)

m_grid  <- if (!is.null(models_BSSSM_2019_beta$m_grid))  models_BSSSM_2019_beta$m_grid  else 100L
bm_grid <- if (!is.null(models_BSSSM_2019_beta$bm_grid)) models_BSSSM_2019_beta$bm_grid else 4
dt      <- if (!is.null(models_BSSSM_2019_beta$dt))      models_BSSSM_2019_beta$dt      else (1/252)

m  <- as.integer(m_grid)
bm <- as.numeric(bm_grid)

# ---- Reconstruct grid, Γ, δ, and emission parameters ----
b_edges <- seq(-bm, bm, length.out = m + 1)
h_grid  <- b_edges[2] - b_edges[1]
bstar   <- 0.5 * (b_edges[-1] + b_edges[-(m + 1)])   # m midpoints

# Γ_ij ≈ P(C_t in bin j | C_{t-1} in bin i)
Gamma <- matrix(0, nrow = m, ncol = m)
for (i in 1:m) {
  Gamma[i, ] <- h_grid * dnorm(bstar, mean = phi_hat * bstar[i], sd = sigZ_hat)
}
# Normalize rows to sum to 1 (numerical stability under truncation)
Gamma <- sweep(Gamma, 1, rowSums(Gamma), "/")

# Stationary distribution for AR(1): C_t ~ N(0, sigZ^2/(1-phi^2))
sd0   <- sigZ_hat / sqrt(1 - phi_hat^2)
delta <- h_grid * dnorm(bstar, mean = 0, sd = sd0)
delta <- delta / sum(delta)

# Emissions
mu_vec    <- mu0_hat + beta_mu_hat * bstar
sigma_vec <- sigma_hat * exp(beta_sigma_hat * bstar)

mean_vec <- (mu_vec - 0.5 * sigma_vec^2) * dt
sd_vec   <- sigma_vec * sqrt(dt)

# ---- Forward algorithm: φ_t(j) = P(C_t ≈ b_j | X^{(t)}) ----
n <- n_all

allprobs <- matrix(1, nrow = n, ncol = m)   # default 1 for NAs
ind      <- which(!is.na(x_all))

for (j in 1:m) {
  allprobs[ind, j] <- dnorm(
    x_all[ind],
    mean = mean_vec[j],
    sd   = sd_vec[j]
  )
}

phi   <- matrix(NA_real_, nrow = n, ncol = m)
scale <- numeric(n)

foo      <- delta * allprobs[1, ]
scale[1] <- sum(foo)
if (!is.finite(scale[1]) || scale[1] <= 0) scale[1] <- 1e-300
phi[1, ] <- foo / scale[1]

for (t in 2:n) {
  foo      <- as.numeric(phi[t - 1, ] %*% Gamma) * allprobs[t, ]
  scale[t] <- sum(foo)
  if (!is.finite(scale[t]) || scale[t] <= 0) scale[t] <- 1e-300
  phi[t, ] <- foo / scale[t]
}

# ---- Mixture helpers: CDF, quantile, mode, stats for 90% intervals ----
mixture_cdf <- function(z, w, mu, sd) sum(w * pnorm(z, mean = mu, sd = sd))

mixture_quantile <- function(p, w, mu, sd) {
  mu_mix  <- sum(w * mu)
  var_mix <- sum(w * (sd^2 + mu^2)) - mu_mix^2
  sd_mix  <- sqrt(max(var_mix, 1e-12))
  lower   <- mu_mix - 8 * sd_mix
  upper   <- mu_mix + 8 * sd_mix
  uniroot(
    f = function(z) mixture_cdf(z, w, mu, sd) - p,
    interval = c(lower, upper)
  )$root
}

mixture_mode <- function(w, mu, sd) {
  mu_mix  <- sum(w * mu)
  var_mix <- sum(w * (sd^2 + mu^2)) - mu_mix^2
  sd_mix  <- sqrt(max(var_mix, 1e-12))
  lower   <- mu_mix - 8 * sd_mix
  upper   <- mu_mix + 8 * sd_mix
  f_pdf   <- function(z) sum(w * dnorm(z, mean = mu, sd = sd))
  optimize(f_pdf, interval = c(lower, upper), maximum = TRUE)$maximum
}

mixture_stats_90 <- function(w, mu, sd) {
  list(
    mean   = sum(w * mu),
    mode   = mixture_mode(w, mu, sd),
    median = mixture_quantile(0.5,  w, mu, sd),
    lo90   = mixture_quantile(0.05, w, mu, sd),
    hi90   = mixture_quantile(0.95, w, mu, sd)
  )
}

# ---- Γ^h and horizons for forecast table ----
mat_power <- function(M, h) {
  if (h == 0) return(diag(nrow(M)))
  R <- diag(nrow(M))
  for (i in 1:h) R <- R %*% M
  R
}

max_h <- n_all - n_train

candidate_h <- c(1L, 5L, 21L, 63L, 252L, 3L * 252L, max_h)
h_vec <- sort(unique(candidate_h[candidate_h > 0 & candidate_h <= max_h]))

Gamma_powers <- lapply(h_vec, function(h) mat_power(Gamma, h))
Delta <- dt

# ---- Forecast table from origin T0 = n_train ----
T0_idx  <- n_train
T0_date <- dates[T0_idx]

rows <- vector("list", length(h_vec))

for (k in seq_along(h_vec)) {
  h  <- h_vec[k]
  Gh <- Gamma_powers[[k]]
  
  w_h <- as.numeric(phi[T0_idx, ] %*% Gh)
  st  <- mixture_stats_90(w_h, mean_vec, sd_vec)
  
  target_idx  <- T0_idx + h
  target_date <- dates[target_idx]
  target_year <- year(target_date)
  
  rows[[k]] <- tibble(
    year                = target_year,
    horizon             = h * Delta,   # in years
    forecast_mode       = st$mode,
    forecast_median     = st$median,
    forecast_mean       = st$mean,
    nominal_90_interval = sprintf("[%.6f, %.6f]", st$lo90, st$hi90)
  )
}

BSSSM_beta_forecast_table_90 <- bind_rows(rows)

# ---- 1-step-ahead out-of-sample prediction errors (mixture-mean predictor) ----
# ξ_t(1) = φ_t Γ ;  E[X_{t+1}|X^{(t)}] = Σ_j ξ_{t,j}(1) * mean_vec[j]
xhat_all <- rep(NA_real_, n)

for (t in 1:(n - 1)) {
  w_1step <- as.numeric(phi[t, ] %*% Gamma)
  xhat_all[t + 1] <- sum(w_1step * mean_vec)
}

xhat_test_ssm <- xhat_all[idx_test]
err_test_ssm  <- x_test - xhat_test_ssm

MSE_ssm_1  <- mean(err_test_ssm^2, na.rm = TRUE)
RMSE_ssm_1 <- sqrt(MSE_ssm_1)
MAE_ssm_1  <- mean(abs(err_test_ssm), na.rm = TRUE)

BSSSM_beta_error_results_1step <- tibble(
  Metric = c("MSE", "RMSE", "MAE"),
  Value  = c(MSE_ssm_1, RMSE_ssm_1, MAE_ssm_1)
)

# ---- Output ----
print(BSSSM_beta_forecast_table_90, n = Inf)
print(BSSSM_beta_error_results_1step, n = Inf)

# ---- LaTeX table output: forecast table ----
# library(knitr); library(kableExtra)
# BSSSM_beta_forecast_table_90 %>%
#   mutate(
#     across(c(horizon, forecast_mode, forecast_median, forecast_mean), ~ round(.x, 6))
#   ) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "BS-SSM$_{\\beta}$: mixture forecasts from origin $T_0$ (end of training sample) with nominal 90\\% mixture intervals.",
#     label    = "tab:BSSSMbetaForecast90"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- LaTeX table output: 1-step-ahead errors ----
# BSSSM_beta_error_results_1step %>%
#   mutate(Value = signif(Value, 8)) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "Out-of-sample prediction errors for BS-SSM$_{\\beta}$ (1-step-ahead mixture-mean predictor, test set).",
#     label    = "tab:BSSSMbetaErrors1Step"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- Keep everything for later comparison / LaTeX ----
BSSSM_beta_results <- list(
  forecast_table    = BSSSM_beta_forecast_table_90,
  errors_1step      = BSSSM_beta_error_results_1step,
  theta_hat_nat     = theta_hat_nat,
  phi_hat           = phi_hat,
  sigZ_hat          = sigZ_hat,
  mu0_hat           = mu0_hat,
  sigma_hat         = sigma_hat,
  beta_mu_hat       = beta_mu_hat,
  beta_sigma_hat    = beta_sigma_hat,
  m_grid            = m,
  bm_grid           = bm,
  dt                = dt,
  bstar             = bstar,
  mean_vec          = mean_vec,
  sd_vec            = sd_vec,
  Gamma             = Gamma,
  delta             = delta,
  x_test            = x_test,
  xhat_test         = xhat_test_ssm,
  residuals         = err_test_ssm,
  test_dates        = date_test,
  horizons_days     = h_vec
)
