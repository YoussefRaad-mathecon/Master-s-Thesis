####################################################################################################################
####################################################################################################################
#---------------------------------------- Fitting BS ----------------------------------------------------------
####################################################################################################################
##################################################################################################################

### Load library, data and WD
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library("fHMM") ### Model fitting
library(Quandl) ### Data
library(dplyr) ### Data manipulations
library(lubridate) ### Time-variables
library(stats) ### ACF plots
library(matrixcalc) ### Matrix calculations
library("RColorBrewer") ### Colors
library(latex2exp) ### Text for plots
library(matrixStats) ### ColSds
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)
load("sp500_clean.Rdata")


### Parameters
dt <- 1/252
spx <- spx %>% slice_head(n = 23109)     

S <- na.omit(spx$Return)
R <- na.omit(spx$Return)
#----------------------------------------------Log-likelihood--------------------------------------------------#
BS_loglik <- function(params, S) {
  mu <- params[1]
  sigma <- exp(params[2])  # keep positivity via log-param like your Vasicek code
  
  n <- length(S)
  # log returns
  R <- na.omit(spx$Return)
  
  mean_R <- (mu - 0.5 * sigma^2) * dt
  var_R  <- sigma^2 * dt
  
  ll <- sum(dnorm(R, mean = mean_R, sd = sqrt(var_R), log = TRUE))
  return(-ll)  # negative log-likelihood for minimization
}

#-------------------------------------------Fitting---------------------------------------------------------------------#
# Continue with the BS log-likelihood optimization as before
mod <- function(S) {
  init_params <- c(0.06, log(0.2))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)
#------------------------------------------------------------------------------------------------------#
mod <- function(S) {
  init_params <- c(0.02, log(0.6))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)
#------------------------------------------------------------------------------------------------------#
# one-line fit (matches your nlm HMM call style)
mod_BS <- nlm(BS_loglik, c(0.06, log(0.2)), S = S, print.level = 2, iterlim = 10000, hessian = TRUE)

# one-line print (k = 2)
cat(sprintf("mu=%.6f  sigma=%.6f  AIC=%.3f  BIC=%.3f\n",
            mod_BS$estimate[1], exp(mod_BS$estimate[2]),
            2*2 + 2*mod_BS$minimum,
            2*log(length(S)) + 2*mod_BS$minimum))


models_BS_2019 <- list(fit_BS)

save(models_BS_2019, file = "fitted_params_BS_2019.RData")


fit_BS
residuals(fit_BS)
# Pieces for AIC/BIC
k          <- 2                         # parameters: mu and sigma
n_eff      <- length(S)                 # sample size used in the likelihood
logLik_hat <- -fit_BS$minimum              # because BS_loglik returns NEGATIVE log-lik

# Compute + print
AIC  <- 2*k - 2*logLik_hat
BIC  <- k*log(n_eff) - 2*logLik_hat

cat(sprintf("AIC: %.3f\nBIC: %.3f\n", AIC, BIC))


################################################################################################
#--------------------------- Closed-form MLEs for Black–Scholes (no optimization) --------------#
################################################################################################
bs_mle <- function(R, dt) {
  R   <- as.numeric(na.omit(R))
  Tn  <- length(R)
  xbar <- mean(R)
  
  # MLEs from your derivation:
  sigma2_hat <- sum((R - xbar)^2) / (Tn * dt)     # \hat{\sigma}^2
  mu_hat     <- xbar / dt + 0.5 * sigma2_hat      # \hat{\mu}
  sigma_hat  <- sqrt(sigma2_hat)
  
  # Log-likelihood at the MLE:
  var_hat    <- sigma2_hat * dt                   # = \hat{\sigma}^2 \Delta
  logLik_hat <- -(Tn/2) * (log(2*pi*var_hat) + 1) # = -T/2 [ log(2π σ^2 Δ) + 1 ]
  
  # Information criteria (k = 2: mu and sigma)
  k   <- 2L
  AIC <- 2*k - 2*logLik_hat
  BIC <- k*log(Tn) - 2*logLik_hat
  
  list(mu = mu_hat,
       sigma = sigma_hat,
       sigma2 = sigma2_hat,
       T = Tn,
       dt = dt,
       xbar = xbar,
       logLik = logLik_hat,
       AIC = AIC,
       BIC = BIC)
}

# Run closed-form MLEs on your data:
fit_BS_mle <- bs_mle(S, dt)

cat(sprintf("Closed-form MLE  ->  mu=%.8f  sigma=%.8f  logLik=%.3f  AIC=%.3f  BIC=%.3f\n",
            fit_BS_mle$mu, fit_BS_mle$sigma, fit_BS_mle$logLik, fit_BS_mle$AIC, fit_BS_mle$BIC))


################################################################################################
#--------------------------- Residuals for Black–Scholes fit ----------------------------------#
################################################################################################

# Use the same return series as in the likelihood
R  <- as.numeric(na.omit(spx$Return))

# --- choose which fit to use ---------------------------------------------------
mu_hat    <- mod_BS$estimate[1]
sigma_hat <- exp(mod_BS$estimate[2])

# --- compute standardized residuals -------------------------------------------
mean_R_hat <- (mu_hat - 0.5 * sigma_hat^2) * dt
sd_R_hat   <- sigma_hat * sqrt(dt)

res_raw <- R - mean_R_hat
res_std <- res_raw / sd_R_hat   # should be ~ N(0,1) under the model

res_df <- tibble(
  t       = seq_along(res_std),
  res_raw = res_raw,
  res_std = res_std
)

# Keep only finite residuals
res_std_finite <- res_std[is.finite(res_std)]

# (Optional) date alignment if you want time on x-axis
# Adjust column name if your date column is different
date_aligned <- spx$Date[!is.na(spx$Return)]  # change "Date" if needed



## 2) Big 2×2 diagnostic figure -----------------------------------------------

png("BS_ResidualDiagnostics.png", width = 2400, height = 2400, res = 300)

# Slightly larger top margin so titles don't get clipped
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.5, 0))

col_main  <- "#901a1E"
col_ref   <- "steelblue"

### (a) Time series of standardized residuals
if (exists("date_aligned") && length(date_aligned) == length(res_std)) {
  plot(date_aligned, res_std, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Standardized residuals")
} else {
  # Fallback: index on x-axis
  plot(seq_along(res_std), res_std, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Standardized residuals")
}
abline(h = 0, lwd = 1, col = "grey60")
abline(h = c(-3, 3), lty = 2, col = "grey75")
title("Time series", cex.main = 1.4)

### (b) Histogram (trimmed to 1–99% quantiles) + N(0,1) overlay
q_low  <- quantile(res_std_finite, 0.01, na.rm = TRUE)
q_high <- quantile(res_std_finite, 0.99, na.rm = TRUE)
res_trimmed <- res_std_finite[res_std_finite >= q_low & res_std_finite <= q_high]

hist(res_trimmed,
     breaks = 40,
     freq   = FALSE,
     col    = col_main,
     border = "white",
     main   = "Histogram (trimmed 1–99%)",
     xlab   = "Standardized residual",
     ylab   = "Density")
curve(dnorm(x, mean = 0, sd = 1),
      add = TRUE, col = col_ref, lwd = 2)

### (c) QQ-plot vs N(0,1)
qqnorm(res_std_finite,
       main = "QQ-plot",
       col  = col_main,
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles",
       cex.main = 1.4)
qqline(res_std_finite, col = col_ref, lwd = 2)

### (d) ACF of standardized residuals
acf_obj <- acf(res_std_finite, plot = FALSE)
plot(acf_obj,
     main = "",                 # leave empty here
     xlab = "Lag",
     ylab = "Autocorrelation",
     col  = col_main,
     lwd  = 2)
title("ACF", cex.main = 1.4)    # force the ACF title

mtext("BSM Residuals",
      side  = 3,
      outer = TRUE,
      cex   = 1,
      font  = 2)   # 2 = bold
dev.off()



####### SE

fit_BS <- mod_BS   # just to match your other naming

# Hessian on working scale: eta = (mu, log sigma)
H_w <- fit_BS$hessian
theta_hat_tr <- fit_BS$estimate  # (mu, log sigma)

if (!all(dim(H_w) == c(2, 2))) stop("Unexpected Hessian dimension (expected 2x2).")

## --- Back-transform to natural scale ------------------------------------------
theta_hat <- c(
  mu    = theta_hat_tr[1],
  sigma = exp(theta_hat_tr[2])
)

## --- Jacobian d(theta)/d(eta) at the MLE --------------------------------------
# d mu / d mu          = 1
# d sigma / d log(sigma) = sigma
M <- diag(c(1, theta_hat["sigma"]))
rownames(M) <- colnames(M) <- names(theta_hat)

## --- Safe inversion helper ----------------------------------------------------
invert_safe <- function(A) {
  tryCatch(solve(A), error = function(e) {
    message("Hessian not positive-definite; using pseudoinverse via MASS::ginv.")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install.packages('MASS').")
    MASS::ginv(A)
  })
}

## --- Covariance on working scale, then natural scale --------------------------
# nlm Hessian is for -log L, so Var_hat(eta_hat) ≈ H_w^{-1}
Cov_w <- invert_safe(H_w)

# Delta method: Cov(theta_hat) = M %*% Cov(eta_hat) %*% t(M)
Cov_n <- M %*% Cov_w %*% t(M)

## --- Standard errors + 95% Wald CIs on natural scale --------------------------
se <- sqrt(pmax(diag(Cov_n), 0))  # guard against tiny negative diagonals
z  <- 1.96

CI_lower <- theta_hat - z * se
CI_upper <- theta_hat + z * se

se_table_BS <- data.frame(
  Parameter = names(theta_hat),
  Estimate  = as.numeric(theta_hat),
  StdError  = as.numeric(se),
  CI95_L    = as.numeric(CI_lower),
  CI95_U    = as.numeric(CI_upper),
  row.names = NULL,
  check.names = FALSE
)

cat("\nBlack–Scholes MLE with Hessian-based SEs (natural scale):\n")
print(transform(se_table_BS,
                Estimate = signif(Estimate, 6),
                StdError = signif(StdError, 6),
                CI95_L   = signif(CI95_L, 6),
                CI95_U   = signif(CI95_U, 6)),
      row.names = FALSE)


