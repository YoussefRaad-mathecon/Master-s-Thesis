# Black–Scholes (S&P 500 to 2019): Fit, closed-form MLE, residual diagnostics, Hessian SEs

# ---- Packages ----
library(tidyverse)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data ----
load("sp500_clean.Rdata")

dt  <- 1/252
spx <- spx %>% slice_head(n = 23109)

S <- na.omit(spx$Return)
R <- na.omit(spx$Return)

# ---- Log-likelihood ----
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

# ---- Fitting (two initializations + one-line fit) ----
mod <- function(S) {
  init_params <- c(0.06, log(0.2))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)

mod <- function(S) {
  init_params <- c(0.02, log(0.6))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)

# one-line fit (matches your nlm HMM call style)
mod_BS <- nlm(BS_loglik, c(0.06, log(0.2)), S = S, print.level = 2, iterlim = 10000, hessian = TRUE)

# define fit_BS (used below for saving + AIC/BIC)
fit_BS <- mod_BS

# one-line print (k = 2)
cat(sprintf("mu=%.6f  sigma=%.6f  AIC=%.3f  BIC=%.3f\n",
            mod_BS$estimate[1], exp(mod_BS$estimate[2]),
            2*2 + 2*mod_BS$minimum,
            2*log(length(S)) + 2*mod_BS$minimum))

# ---- Save fit ----
models_BS_2019 <- list(fit_BS)
save(models_BS_2019, file = "fitted_params_BS_2019.RData")

fit_BS
residuals(fit_BS)

# ---- AIC/BIC from nlm fit (same likelihood definition) ----
k          <- 2                         # parameters: mu and sigma
n_eff      <- length(S)                 # sample size used in the likelihood
logLik_hat <- -fit_BS$minimum           # because BS_loglik returns NEGATIVE log-lik

AIC  <- 2*k - 2*logLik_hat
BIC  <- k*log(n_eff) - 2*logLik_hat

cat(sprintf("AIC: %.3f\nBIC: %.3f\n", AIC, BIC))


# ---- Closed-form MLEs (no optimization) ----
bs_mle <- function(R, dt) {
  R   <- as.numeric(na.omit(R))
  Tn  <- length(R)
  xbar <- mean(R)
  
  sigma2_hat <- sum((R - xbar)^2) / (Tn * dt)     # \hat{\sigma}^2
  mu_hat     <- xbar / dt + 0.5 * sigma2_hat      # \hat{\mu}
  sigma_hat  <- sqrt(sigma2_hat)
  
  var_hat    <- sigma2_hat * dt
  logLik_hat <- -(Tn/2) * (log(2*pi*var_hat) + 1)
  
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

fit_BS_mle <- bs_mle(S, dt)

cat(sprintf("Closed-form MLE  ->  mu=%.8f  sigma=%.8f  logLik=%.3f  AIC=%.3f  BIC=%.3f\n",
            fit_BS_mle$mu, fit_BS_mle$sigma, fit_BS_mle$logLik, fit_BS_mle$AIC, fit_BS_mle$BIC))


# ---- Residuals (standardized) + diagnostics figure ----
R  <- as.numeric(na.omit(spx$Return))

mu_hat    <- mod_BS$estimate[1]
sigma_hat <- exp(mod_BS$estimate[2])

mean_R_hat <- (mu_hat - 0.5 * sigma_hat^2) * dt
sd_R_hat   <- sigma_hat * sqrt(dt)

res_raw <- R - mean_R_hat
res_std <- res_raw / sd_R_hat

res_df <- tibble(
  t       = seq_along(res_std),
  res_raw = res_raw,
  res_std = res_std
)

res_std_finite <- res_std[is.finite(res_std)]

date_aligned <- spx$Date[!is.na(spx$Return)]  # change "Date" if needed


png("BS_ResidualDiagnostics.png", width = 2400, height = 2400, res = 300)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.5, 0))

col_main  <- "#901a1E"
col_ref   <- "steelblue"

if (exists("date_aligned") && length(date_aligned) == length(res_std)) {
  plot(date_aligned, res_std, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Standardized residuals")
} else {
  plot(seq_along(res_std), res_std, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Standardized residuals")
}
abline(h = 0, lwd = 1, col = "grey60")
abline(h = c(-3, 3), lty = 2, col = "grey75")
title("Time series", cex.main = 1.4)

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

qqnorm(res_std_finite,
       main = "QQ-plot",
       col  = col_main,
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles",
       cex.main = 1.4)
qqline(res_std_finite, col = col_ref, lwd = 2)

acf_obj <- acf(res_std_finite, plot = FALSE)
plot(acf_obj,
     main = "",
     xlab = "Lag",
     ylab = "Autocorrelation",
     col  = col_main,
     lwd  = 2)
title("ACF", cex.main = 1.4)

mtext("BSM Residuals",
      side  = 3,
      outer = TRUE,
      cex   = 1,
      font  = 2)
dev.off()


# ---- Hessian-based SEs (delta method back to natural scale) ----
fit_BS <- mod_BS

H_w <- fit_BS$hessian
theta_hat_tr <- fit_BS$estimate  # (mu, log sigma)

if (!all(dim(H_w) == c(2, 2))) stop("Unexpected Hessian dimension (expected 2x2).")

theta_hat <- c(
  mu    = theta_hat_tr[1],
  sigma = exp(theta_hat_tr[2])
)

M <- diag(c(1, theta_hat["sigma"]))
rownames(M) <- colnames(M) <- names(theta_hat)

invert_safe <- function(A) {
  tryCatch(solve(A), error = function(e) {
    message("Hessian not positive-definite; using pseudoinverse via MASS::ginv.")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install.packages('MASS').")
    MASS::ginv(A)
  })
}

Cov_w <- invert_safe(H_w)
Cov_n <- M %*% Cov_w %*% t(M)

se <- sqrt(pmax(diag(Cov_n), 0))
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
