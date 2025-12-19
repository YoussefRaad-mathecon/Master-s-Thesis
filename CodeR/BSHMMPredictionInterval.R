# 4-state BS–HMM (μ_i, σ_i): Forecast table + 1-step-ahead MSE/RMSE/MAE

# ---- Packages ----
library(tidyverse)
library(lubridate)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data, model, and split ----
load("sp500_clean.Rdata")               # expects spx$Date, spx$Return
load("fitted_params_BSHMM_2019.RData")  # expects models_BSHMM_2019$mod4_mu_sigma

if ("Date" %in% names(spx)) {
  spx <- spx %>% arrange(Date)
}

Delta <- 1/252

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

# ---- Unpack BS–HMM parameters (μ_i, σ_i, Γ, δ) ----
unpack_bshmm <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
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
  } else {  # "mu_sigma"
    mu        <- th[1:N]
    log_sigma <- th[(N + 1):(2 * N)]
    eta       <- th[(2 * N + 1):length(th)]
    sigma     <- exp(log_sigma)
  }
  
  # Transition matrix Γ
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- sweep(Gamma, 1, rowSums(Gamma), "/")
  
  # Stationary distribution δ
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  delta <- as.numeric(delta / sum(delta))
  
  list(mu = mu, sigma = sigma, Gamma = Gamma, delta = delta)
}

if (is.null(models_BSHMM_2019$mod4_mu_sigma)) {
  stop("models_BSHMM_2019$mod4_mu_sigma not found in fitted_params_BSHMM_2019.RData.")
}

N             <- 4L
mod4_mu_sigma <- models_BSHMM_2019$mod4_mu_sigma

pars4 <- unpack_bshmm(mod4_mu_sigma, N = N, family = "mu_sigma")

mu_hat_i    <- pars4$mu
sigma_hat_i <- pars4$sigma
Gamma_hat   <- pars4$Gamma
delta_hat   <- pars4$delta

# State-wise daily mean/sd of log-returns:
# X_t | S_t = i ~ N( (μ_i - 0.5 σ_i^2)*Δ , σ_i^2 * Δ )
mu_state_daily    <- (mu_hat_i - 0.5 * sigma_hat_i^2) * Delta
sigma_state_daily <- sigma_hat_i * sqrt(Delta)

# ---- Forward algorithm: φ_t = P(S_t = i | X^{(t)}) ----
n <- n_all

allprobs <- matrix(1, nrow = n, ncol = N)  # default 1 for NAs
ind      <- which(!is.na(x_all))

for (j in 1:N) {
  allprobs[ind, j] <- dnorm(
    x_all[ind],
    mean = mu_state_daily[j],
    sd   = sigma_state_daily[j]
  )
}

phi   <- matrix(NA_real_, nrow = n, ncol = N)
scale <- numeric(n)

# t = 1
foo      <- delta_hat * allprobs[1, ]
scale[1] <- sum(foo)
phi[1, ] <- foo / scale[1]

# t = 2,...,n
for (t in 2:n) {
  foo      <- (phi[t - 1, ] %*% Gamma_hat) * allprobs[t, ]
  scale[t] <- sum(foo)
  phi[t, ] <- foo / scale[t]
}

# ---- Mixture helpers: CDF, quantile, mode, stats for 90% intervals ----
mixture_cdf <- function(z, w, mu, sd) {
  sum(w * pnorm(z, mean = mu, sd = sd))
}

mixture_quantile <- function(p, w, mu, sd) {
  mu_mix  <- sum(w * mu)
  var_mix <- sum(w * (sd^2 + mu^2)) - mu_mix^2
  sd_mix  <- sqrt(max(var_mix, 1e-12))
  
  lower <- mu_mix - 8 * sd_mix
  upper <- mu_mix + 8 * sd_mix
  
  uniroot(
    f = function(z) mixture_cdf(z, w, mu, sd) - p,
    interval = c(lower, upper)
  )$root
}

mixture_mode <- function(w, mu, sd) {
  mu_mix  <- sum(w * mu)
  var_mix <- sum(w * (sd^2 + mu^2)) - mu_mix^2
  sd_mix  <- sqrt(max(var_mix, 1e-12))
  
  lower <- mu_mix - 8 * sd_mix
  upper <- mu_mix + 8 * sd_mix
  
  f_pdf <- function(z) sum(w * dnorm(z, mean = mu, sd = sd))
  optimize(f_pdf, interval = c(lower, upper), maximum = TRUE)$maximum
}

mixture_stats_90 <- function(w, mu, sd) {
  mu_mix <- sum(w * mu)
  list(
    mean   = mu_mix,
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

max_h <- n_all - n_train   # length of test window

# Same base horizons as BSM, plus full horizon
candidate_h <- c(1L, 5L, 21L, 63L, 252L, 3L * 252L, max_h)
h_vec <- sort(unique(candidate_h[candidate_h > 0 & candidate_h <= max_h]))

Gamma_powers <- lapply(h_vec, function(h) mat_power(Gamma_hat, h))

# ---- Forecast table from origin T0 = n_train ----
T0_idx  <- n_train
T0_date <- dates[T0_idx]

rows <- vector("list", length(h_vec))

for (k in seq_along(h_vec)) {
  h  <- h_vec[k]
  Gh <- Gamma_powers[[k]]
  
  # Forecast weights at horizon h from origin T0: ξ_T0(h) = φ_T0 Γ^h
  w_h <- as.numeric(phi[T0_idx, ] %*% Gh)
  
  st <- mixture_stats_90(w_h, mu_state_daily, sigma_state_daily)
  
  target_idx  <- T0_idx + h
  target_date <- dates[target_idx]
  target_year <- year(target_date)
  
  rows[[k]] <- tibble(
    year                = target_year,
    horizon             = h * Delta,        # in years
    forecast_mode       = st$mode,
    forecast_median     = st$median,
    forecast_mean       = st$mean,
    nominal_90_interval = sprintf("[%.6f, %.6f]", st$lo90, st$hi90)
  )
}

BSHMM_forecast_table_90 <- bind_rows(rows)

# ---- 1-step-ahead out-of-sample prediction errors (mixture-mean predictor) ----
# For each t = 1,...,n-1:
#   ξ_t(1) = φ_t Γ_hat
#   E[X_{t+1} | X^{(t)}] = Σ_i ξ_{t,i}(1) * μ_i^(X),
# where μ_i^(X) = mu_state_daily[i].
xhat_all <- rep(NA_real_, n)  # forecast for X_t from info up to t-1

for (t in 1:(n - 1)) {
  w_1step <- as.numeric(phi[t, ] %*% Gamma_hat)    # ξ_t(1)
  xhat_all[t + 1] <- sum(w_1step * mu_state_daily) # mixture mean
}

xhat_test_hmm <- xhat_all[idx_test]
err_test_hmm  <- x_test - xhat_test_hmm

MSE_hmm_1  <- mean(err_test_hmm^2, na.rm = TRUE)
RMSE_hmm_1 <- sqrt(MSE_hmm_1)
MAE_hmm_1  <- mean(abs(err_test_hmm), na.rm = TRUE)

BSHMM_error_results_1step <- tibble(
  Metric = c("MSE", "RMSE", "MAE"),
  Value  = c(MSE_hmm_1, RMSE_hmm_1, MAE_hmm_1)
)

# ---- Output ----
print(BSHMM_forecast_table_90, n = Inf)
print(BSHMM_error_results_1step, n = Inf)

# ---- LaTeX table output: forecast table ----
# library(knitr); library(kableExtra)
# BSHMM_forecast_table_90 %>%
#   mutate(
#     across(c(horizon, forecast_mode, forecast_median, forecast_mean), ~ round(.x, 6))
#   ) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "4-state BS--HMM $(\\mu_i,\\,\\sigma_i)$: mixture forecasts from origin $T_0$ (end of training sample) with nominal 90\\% mixture intervals.",
#     label    = "tab:BSHMMForecast90"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- LaTeX table output: 1-step-ahead errors ----
# BSHMM_error_results_1step %>%
#   mutate(Value = signif(Value, 8)) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "Out-of-sample prediction errors for the 4-state BS--HMM (1-step-ahead mixture-mean predictor, test set).",
#     label    = "tab:BSHMMErrors1Step"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- Keep everything ----
BSHMM_results <- list(
  forecast_table     = BSHMM_forecast_table_90,
  errors_1step       = BSHMM_error_results_1step,
  mu_hat_i           = mu_hat_i,
  sigma_hat_i        = sigma_hat_i,
  Gamma_hat          = Gamma_hat,
  delta_hat          = delta_hat,
  mu_state_daily     = mu_state_daily,
  sigma_state_daily  = sigma_state_daily,
  x_test             = x_test,
  xhat_test          = xhat_test_hmm,
  residuals          = err_test_hmm,
  test_dates         = date_test
)
