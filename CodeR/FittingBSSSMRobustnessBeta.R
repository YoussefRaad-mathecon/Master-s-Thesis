# BS–SSM_beta (S&P 500 to 2019): Robustness grid over (m, bm)

# ---- Packages ----
library(tidyverse)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)
load("sp500_clean.Rdata")

# ---- Data ----
spx   <- spx %>% slice_head(n = 23109)
S     <- na.omit(spx$Return)    
x_raw <- S
x     <- as.numeric(x_raw)       

# ---- Parameters ----
dt <- 1/252

# ---------------------------------------------------------------------------------------------#
# Likelihood for BS-SSM _beta (one AR(1) state driving both drift and volatility)
#
# State:
#   C_t = phi_state * C_{t-1} + eta_t,  eta_t ~ N(0, sigma_state^2)
#
# Emission:
#   mu_t    = mu0 + beta_mu * C_t
#   sigma_t = sigma * exp(beta_sigma * C_t)
#   X_t|C_t ~ N( (mu_t - 0.5 * sigma_t^2) dt,  sigma_t^2 dt )
#
# Parameterization of theta.star:
#   theta.star[1] = logit(phi_state)     -> phi_state in (0,1)
#   theta.star[2] = log(sigma_state)     -> sigma_state > 0
#   theta.star[3] = mu0                  -> R
#   theta.star[4] = log(sigma)           -> sigma > 0
#   theta.star[5] = beta_mu              -> R
#   theta.star[6] = beta_sigma           -> R
# ---------------------------------------------------------------------------------------------#
mllk_BS_beta <- function(theta.star, x, m, bm, dt = 1/252) {
  # Transformed parameters
  phi_state   <- plogis(theta.star[1])   # phi in (0,1)
  sigma_state <- exp(theta.star[2])      # >0
  mu0         <- theta.star[3]           # baseline drift
  sigma       <- exp(theta.star[4])      # baseline vol >0
  beta_mu     <- theta.star[5]           # loading in drift
  beta_sigma  <- theta.star[6]           # loading in log-vol
  
  # Latent state grid for C_t
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])      # midpoints
  
  # Transition kernel for C_t (discretized Gaussian AR(1))
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(
      bstar,
      mean = phi_state * bstar[i],
      sd   = sigma_state
    )
  }
  
  # Stationary initial distribution for C_t
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  
  # Emission (returns conditional on C_t ≈ bstar)
  mu_vec      <- mu0 + beta_mu * bstar
  sigma_vec   <- sigma * exp(beta_sigma * bstar)
  
  mean_vec    <- (mu_vec - 0.5 * sigma_vec^2) * dt
  sd_vec      <- sigma_vec * sqrt(dt)
  
  # Forward algorithm
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
  
  return(-logL)  # nlm minimizes
}

# ---- Initial guess ----
theta.start_beta <- c(
  qlogis(0.95),  # phi_state
  log(0.08),     # sigma_state
  0.08,          # mu0 (baseline annual drift ~8%)
  log(0.15),     # sigma (baseline vol ~15% p.a.)
  0.00,          # beta_mu
  0.50           # beta_sigma
)

# ---- Robustness grid ----
m_values  <- c(20, 40, 70, 100, 200)
bm_values <- c(0.5, 1, 2, 3, 4)
combo_grid <- expand.grid(m = m_values, bm = bm_values, KEEP.OUT.ATTRS = FALSE)

# Containers
models_BSSSM_beta_grid <- vector("list", nrow(combo_grid))
names(models_BSSSM_beta_grid) <- sprintf("m%g_bm%g", combo_grid$m, combo_grid$bm)

grid_summary_beta <- tibble::tibble(
  name       = names(models_BSSSM_beta_grid),
  m          = combo_grid$m,
  bm         = combo_grid$bm,
  status     = NA_character_,
  conv_code  = NA_integer_,
  negloglik  = NA_real_,
  phi        = NA_real_,
  sigZ       = NA_real_,
  mu0        = NA_real_,
  sigma      = NA_real_,
  beta_mu    = NA_real_,
  beta_sigma = NA_real_,
  elapsed_s  = NA_real_,
  message    = NA_character_
)

# Back-transform helper for the _beta model
bt_beta <- function(est) c(
  phi        = plogis(est[1]),
  sigZ       = exp(est[2]),
  mu0        = est[3],
  sigma      = exp(est[4]),
  beta_mu    = est[5],
  beta_sigma = est[6]
)

# Loop over all (m, bm) combos
for (i in seq_len(nrow(combo_grid))) {
  mi  <- combo_grid$m[i]
  bmi <- combo_grid$bm[i]
  t0  <- Sys.time()
  cat(sprintf("Fitting _beta combo %2d/%2d: m = %-3d, bm = %-3.1f ... ",
              i, nrow(combo_grid), mi, bmi))
  
  fit_i <- try(
    nlm(
      f = mllk_BS_beta,
      p = theta.start_beta,
      x = x,
      m = mi,
      bm = bmi,
      dt = dt,
      print.level = 0,
      hessian = TRUE
    ),
    silent = TRUE
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  if (inherits(fit_i, "try-error")) {
    msg <- as.character(fit_i)
    grid_summary_beta$status[i]    <- "fail"
    grid_summary_beta$message[i]   <- msg
    grid_summary_beta$elapsed_s[i] <- elapsed
    
    models_BSSSM_beta_grid[[i]] <- list(
      status     = "fail",
      m          = mi,
      bm         = bmi,
      error      = msg,
      elapsed_s  = elapsed
    )
    cat("FAIL\n")
    next
  }
  
  theta_hat_tr_i <- fit_i$estimate
  theta_hat_i    <- bt_beta(theta_hat_tr_i)
  
  grid_summary_beta$status[i]     <- "ok"
  grid_summary_beta$conv_code[i]  <- fit_i$code
  grid_summary_beta$negloglik[i]  <- fit_i$minimum
  grid_summary_beta$phi[i]        <- theta_hat_i["phi"]
  grid_summary_beta$sigZ[i]       <- theta_hat_i["sigZ"]
  grid_summary_beta$mu0[i]        <- theta_hat_i["mu0"]
  grid_summary_beta$sigma[i]      <- theta_hat_i["sigma"]
  grid_summary_beta$beta_mu[i]    <- theta_hat_i["beta_mu"]
  grid_summary_beta$beta_sigma[i] <- theta_hat_i["beta_sigma"]
  grid_summary_beta$elapsed_s[i]  <- elapsed
  
  models_BSSSM_beta_grid[[i]] <- list(
    status        = "ok",
    m             = mi,
    bm            = bmi,
    fit           = fit_i,
    theta_tr      = theta_hat_tr_i,
    theta_bt      = theta_hat_i,
    negloglik     = fit_i$minimum,
    conv_code     = fit_i$code,
    hessian       = fit_i$hessian,
    elapsed_s     = elapsed
  )
  
  cat(sprintf("OK (code=%d, -logL=%.6f, %.2fs)\n",
              fit_i$code, fit_i$minimum, elapsed))
}

# ---- Save robustness results ----
models_BSSSM_robustness_beta <- list(
  grid    = models_BSSSM_beta_grid,
  summary = grid_summary_beta
)

grid_summary_beta_sorted <- dplyr::arrange(grid_summary_beta, bm, m)
print(grid_summary_beta_sorted, n = Inf)

save(
  models_BSSSM_robustness_beta,
  grid_summary_beta,
  grid_summary_beta_sorted,
  file = "fitted_params_BSSSM_grid_2019_beta.RData"
)

saveRDS(grid_summary_beta,        "BSSSM_grid_summary_2019_beta.rds")
saveRDS(grid_summary_beta_sorted, "BSSSM_grid_summary_sorted_2019_beta.rds")

# ---- Load and display saved results ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("fitted_params_BSSSM_grid_2019_beta.RData"); print(grid_summary_beta_sorted, n = Inf)
print(readRDS("BSSSM_grid_summary_sorted_2019_beta.rds"), n = Inf)

options(digits = 22, scipen = 999, pillar.sigfig = 17, tibble.width = Inf)

fmt <- function(x, k = 16) formatC(x, format = "f", digits = k, big.mark = ",", decimal.mark = ".")
grid_summary_beta_sorted %>%
  arrange(negloglik) %>%
  mutate(negloglik_full = fmt(negloglik)) %>%
  select(name, m, bm, conv_code, negloglik_full) %>%
  print(n = Inf, width = Inf)
