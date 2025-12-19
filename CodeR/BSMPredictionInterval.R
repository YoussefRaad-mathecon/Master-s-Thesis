# Constant BSM: Forecast table + 1-step-ahead MSE/RMSE/MAE

# ---- Packages ----
library(tidyverse)
library(lubridate)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data and train/test split ----
load("sp500_clean.Rdata")   # expects spx$Date, spx$Return

# Ensure sorted by date
if ("Date" %in% names(spx)) {
  spx <- spx %>% arrange(Date)
}

x_all <- as.numeric(spx$Return)
dates <- if ("Date" %in% names(spx)) spx$Date else seq_along(x_all)
n_all <- length(x_all)

# Training length as in your other scripts
n_train <- 23109L
if (n_all <= n_train) stop("Not enough observations for a test set: n_total <= 23109.")

idx_train <- 1:n_train
idx_test  <- (n_train + 1):n_all

x_train <- x_all[idx_train]
x_test  <- x_all[idx_test]

date_train <- dates[idx_train]
date_test  <- dates[idx_test]

# ---- Fit iid Normal model on training returns (BSM on log-returns) ----
m_hat  <- mean(x_train, na.rm = TRUE)   # daily mean
v_hat  <- var(x_train,  na.rm = TRUE)   # daily variance
sd_hat <- sqrt(v_hat)

# For a Normal: mode = median = mean = m_hat
forecast_mode   <- m_hat
forecast_median <- m_hat
forecast_mean   <- m_hat

# Nominal 90% interval for a single daily return
z_90     <- qnorm(0.95)
lower_90 <- m_hat - z_90 * sd_hat
upper_90 <- m_hat + z_90 * sd_hat

# ---- Zucchini-style forecast table for horizons up to full period ----
Delta <- 1/252   # trading-year conversion

max_h <- n_all - n_train   # length of test window in trading days

candidate_h <- c(1L, 5L, 21L, 63L, 252L, 3L * 252L, max_h)
h_vec <- sort(unique(candidate_h[candidate_h > 0 & candidate_h <= max_h]))

# Origin for table: last training day
T0_idx  <- n_train
T0_date <- dates[T0_idx]

rows <- vector("list", length(h_vec))

for (k in seq_along(h_vec)) {
  h <- h_vec[k]
  
  target_idx  <- T0_idx + h
  target_date <- dates[target_idx]
  target_year <- year(target_date)
  
  rows[[k]] <- tibble(
    year                = target_year,
    horizon             = h * Delta,   # in years
    forecast_mode       = forecast_mode,
    forecast_median     = forecast_median,
    forecast_mean       = forecast_mean,
    nominal_90_interval = sprintf("[%.6f, %.6f]", lower_90, upper_90)
  )
}

BSM_forecast_table_90 <- bind_rows(rows)

# ---- 1-step-ahead out-of-sample prediction errors (test set) ----
# Under iid Normal BSM, the optimal 1-step-ahead predictor for a daily return is m_hat.
xhat_test_bsm <- rep(m_hat, length(x_test))
err_test_bsm  <- x_test - xhat_test_bsm

MSE_bsm  <- mean(err_test_bsm^2, na.rm = TRUE)
RMSE_bsm <- sqrt(MSE_bsm)
MAE_bsm  <- mean(abs(err_test_bsm), na.rm = TRUE)

BSM_error_results <- tibble(
  Metric = c("MSE", "RMSE", "MAE"),
  Value  = c(MSE_bsm, RMSE_bsm, MAE_bsm)
)

# ---- Output ----
print(BSM_forecast_table_90, n = Inf)
print(BSM_error_results, n = Inf)

# ---- LaTeX table output: forecast table ----
# library(knitr); library(kableExtra)
# BSM_forecast_table_90 %>%
#   mutate(
#     across(c(horizon, forecast_mode, forecast_median, forecast_mean), ~ round(.x, 6))
#   ) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "Constant BSM: forecasts from origin $T_0$ (end of training sample) with nominal 90\\% intervals under an iid Normal model for daily log-returns.",
#     label    = "tab:BSMForecast90"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- LaTeX table output: 1-step-ahead errors ----
# BSM_error_results %>%
#   mutate(Value = signif(Value, 8)) %>%
#   kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "Out-of-sample prediction errors for the constant BSM (1-step-ahead, test set).",
#     label    = "tab:BSMErrors1Step"
#   ) %>%
#   kable_styling(position = "center", latex_options = "hold_position")

# ---- Keep objects for later comparison / LaTeX ----
BSM_results <- list(
  forecast_table = BSM_forecast_table_90,
  errors         = BSM_error_results,
  m_hat          = m_hat,
  sd_hat         = sd_hat,
  date_test      = date_test,
  x_test         = x_test,
  xhat_test      = xhat_test_bsm,
  residuals      = err_test_bsm
)
