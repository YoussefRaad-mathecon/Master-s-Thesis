# BS–SSM_beta (S&P 500 to 2019): Fit + Hessian-based SE diagnostics

# ---- Packages ----
library(tidyverse)
library(MASS)   # ginv() fallback if Hessian is not invertible

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data ----
load("sp500_clean.Rdata")

spx <- spx %>% slice_head(n = 23109)

S     <- na.omit(spx$Return)
x_raw <- S
x     <- as.numeric(x_raw)

summary(x)

# ---- Controls ----
dt      <- 1/252
m_grid  <- 100
bm_grid <- 4

# ---- Likelihood: BS–SSM_beta ----
mllk_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  phi_state    <- plogis(theta.star[1])
  sigma_state  <- exp(theta.star[2])
  mu0          <- theta.star[3]
  sigma        <- exp(theta.star[4])
  beta_mu      <- theta.star[5]
  beta_sigma   <- theta.star[6]
  
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  
  mu_vec      <- mu0 + beta_mu * bstar
  exp_term    <- exp(beta_sigma * bstar)
  sigma_vec   <- sigma * exp_term
  
  mean_vec    <- (mu_vec - 0.5 * sigma_vec^2) * dt
  sd_vec      <- sigma_vec * sqrt(dt)
  
  dens1 <- dnorm(x[1], mean = mean_vec, sd = sd_vec)
  foo   <- delta * dens1
  s1    <- sum(foo)
  logL  <- log(s1)
  phi   <- foo / s1
  
  for (t in 2:length(x)) {
    denst <- dnorm(x[t], mean = mean_vec, sd = sd_vec)
    foo   <- as.numeric(phi %*% Gamma) * denst
    st    <- sum(foo)
    logL  <- logL + log(st)
    phi   <- foo / st
  }
  
  return(-logL)
}

# ---- Initial guess ----
theta.start <- c(
  qlogis(0.95),
  log(0.08),
  0.08,
  log(0.15),
  0.00,
  0.50
)

# ---- Fit ----
fit_BS <- nlm(
  f           = mllk_BS,
  p           = theta.start,
  x           = x,
  m           = m_grid,
  bm          = bm_grid,
  dt          = dt,
  print.level = 2,
  hessian     = TRUE
)

# ---- Back-transform ----
theta_hat_tr <- fit_BS$estimate

theta_BS_hat <- c(
  phi        = plogis(theta_hat_tr[1]),
  sigZ       = exp(theta_hat_tr[2]),
  mu0        = theta_hat_tr[3],
  sigma      = exp(theta_hat_tr[4]),
  beta_mu    = theta_hat_tr[5],
  beta_sigma = theta_hat_tr[6]
)

cat("\n--- Natural-scale parameter estimates (theta_BS_hat) ---\n")
print(theta_BS_hat)

# ---- Hessian-based SE diagnostics ----
H_w <- fit_BS$hessian
cat("\nDiagonal of working-scale Hessian H_w:\n")
print(diag(H_w))

invert_safe_noclip <- function(A) {
  tryCatch(
    solve(A),
    error = function(e) {
      message("solve(H) failed; using MASS::ginv. Error: ", conditionMessage(e))
      MASS::ginv(A)
    }
  )
}

Cov_w <- invert_safe_noclip(H_w)
cat("\nDiagonal of Cov_w (working-scale covariance):\n")
print(diag(Cov_w))

phi_hat        <- theta_BS_hat["phi"]
sigZ_hat       <- theta_BS_hat["sigZ"]
mu0_hat        <- theta_BS_hat["mu0"]
sigma_hat      <- theta_BS_hat["sigma"]
beta_mu_hat    <- theta_BS_hat["beta_mu"]
beta_sigma_hat <- theta_BS_hat["beta_sigma"]

theta_natural <- c(
  phi        = phi_hat,
  sigZ       = sigZ_hat,
  mu0        = mu0_hat,
  sigma      = sigma_hat,
  beta_mu    = beta_mu_hat,
  beta_sigma = beta_sigma_hat
)

J_diag <- c(
  phi_hat   * (1 - phi_hat),
  sigZ_hat,
  1,
  sigma_hat,
  1,
  1
)

M <- diag(J_diag)
rownames(M) <- colnames(M) <- names(theta_natural)

cat("\nJacobian diagonal (natural wrt working):\n")
print(J_diag)

Cov_n_raw <- M %*% Cov_w %*% t(M)
var_n_raw <- diag(Cov_n_raw)
se_raw    <- sqrt(var_n_raw)

cat("\nDiagonal of Cov_n_raw (natural-scale variances, before pmax):\n")
print(var_n_raw)

cat("\nRaw SEs on natural scale (may be NaN if variance < 0):\n")
print(se_raw)

var_clipped <- pmax(var_n_raw, 0)
se_clipped  <- sqrt(var_clipped)

cat("\nClipped variances (pmax(var_n_raw, 0)):\n")
print(var_clipped)

cat("\nClipped SEs (sqrt of clipped variances):\n")
print(se_clipped)

se_summary <- data.frame(
  parameter   = names(theta_natural),
  estimate    = as.numeric(theta_natural),
  var_raw     = as.numeric(var_n_raw),
  se_raw      = as.numeric(se_raw),
  var_clipped = as.numeric(var_clipped),
  se_clipped  = as.numeric(se_clipped),
  row.names   = NULL
)

cat("\n--- Hessian-based SE summary for BS-SSM_beta (natural scale) ---\n")
print(se_summary)

cat("\n--- Focus on beta_mu and beta_sigma ---\n")
print(se_summary[se_summary$parameter %in% c("beta_mu", "beta_sigma"), ])

eigs <- eigen(H_w, symmetric = TRUE, only.values = TRUE)$values
cat("\nEigenvalues of H_w (working-scale Hessian):\n")
print(eigs)
cat("\nCondition ratio (max/min abs eigenvalue):\n")
cat(max(abs(eigs)) / min(abs(eigs)), "\n")

# ---- Save ----
models_BSSSM_2019_beta <- list(
  fit_object  = fit_BS,
  theta_hat   = theta_BS_hat,
  se_summary  = se_summary,
  Cov_w       = Cov_w,
  Cov_n_raw   = Cov_n_raw
)

save(models_BSSSM_2019_beta, file = "fitted_params_BSSSM_2019_beta.RData")
cat("\nSaved models_BSSSM_2019_beta to fitted_params_BSSSM_2019_beta.RData\n")
