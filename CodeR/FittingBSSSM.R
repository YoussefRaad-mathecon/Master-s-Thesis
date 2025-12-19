# BS–SSM (S&P 500 to 2019): Fit Black–Scholes SSM (GBM) on a discretized AR(1) latent grid

# ---- Packages ----
library(tidyverse)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data ----
load("sp500_clean.Rdata")
spx <- spx %>% slice_head(n = 23109)
S <- na.omit(spx$Return)
x_raw <- S           # already log returns
x     <- as.numeric(x_raw)   # coerce to numeric vector

# ---- Controls ----
dt      <- 1/252
m_grid  <- 70
bm_grid <- 3

# ---- Likelihood: BS-SSM ----
mllk_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  phi_state   <- plogis(theta.star[1])      # (0,1)
  sigma_state <- exp(theta.star[2])         # >0
  mu          <- theta.star[3]              # real
  sigma       <- exp(theta.star[4])         # >0
  
  # Latent state grid
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])
  
  # Transition kernel for Z_t
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  # Stationary init
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  
  # Emission
  expZ     <- exp(bstar)
  mean_vec <- (mu * expZ - 0.5 * (sigma * expZ)^2) * dt
  sd_vec   <- (sigma * expZ) * sqrt(dt)
  
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
  
  return(-logL)
}

# ---- Initial guess ----
theta.start <- c(
  qlogis(0.95),  # phi
  log(0.08),     # sigma_state
  0.3,          # mu
  log(0.08)      # sigma
)

# ---- Optimize ----
fit_BS <- nlm(
  f = mllk_BS,
  p = theta.start,
  x = x,
  m = m_grid,
  bm = bm_grid,
  dt = dt,
  print.level = 2, hessian = TRUE
)

# ---- Back-transform to natural scale ----
theta_hat_tr <- fit_BS$estimate
theta_BS_hat <- c(
  phi   = plogis(theta_hat_tr[1]),
  sigZ  = exp(theta_hat_tr[2]),
  mu    = theta_hat_tr[3],
  sigma = exp(theta_hat_tr[4])
)
theta_BS_hat

# ---- Save fit ----
models_BSSSM_2019 <- list(fit_BS)

save(models_BSSSM_2019, file = "fitted_params_BSSSM_2019.RData")
