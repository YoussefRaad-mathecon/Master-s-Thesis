####################################################################################################################
####################################################################################################################
#----------------------------- Continuous State Space Simulation: Black–Scholes (GBM) ------------------------------
####################################################################################################################
####################################################################################################################

## Time step (use 1/252 for daily data)
dt <- 1/252

## True parameters
n             <- 10000
mu_true       <- 0.06          # annualized drift
sigma_true    <- 0.02          # annualized vol
phi_s_true    <- 0.98          # AR(1) persistence of latent state
sigma_s_true  <- 0.02          # AR(1) innovation s.d. of latent state
S0            <- 100

##------------------------------------------------------------------------------
## Simulation of GBM with latent AR(1) state that scales drift and vol:
##   Z_t = phi * Z_{t-1} + sigma_s * eps_t
##   r_t = log S_t - log S_{t-1} ~ N( (mu*e^{Z_t} - 0.5 (sigma*e^{Z_t})^2) dt, (sigma*e^{Z_t})^2 dt )
##------------------------------------------------------------------------------
simulate_SSM_BS <- function(n, mu_true, sigma_true, phi_s_true, sigma_s_true, S0, dt = 1/252) {
  Z      <- numeric(n)            # latent state
  logS   <- numeric(n)            # log prices
  r      <- numeric(n - 1)        # log-returns
  
  Z[1]    <- 0.0
  logS[1] <- log(S0)
  
  for (i in 2:n) {
    # AR(1) latent state
    Z[i] <- phi_s_true * Z[i - 1] + sigma_s_true * rnorm(1)
    
    # State-dependent drift and vol
    mu_i    <- mu_true    * exp(Z[i])
    sigma_i <- sigma_true * exp(Z[i])
    
    # GBM log increment
    incr   <- (mu_i - 0.5 * sigma_i^2) * dt + sigma_i * sqrt(dt) * rnorm(1)
    logS[i] <- logS[i - 1] + incr
    r[i - 1] <- incr
  }
  
  S <- exp(logS)
  list(r = r, Z = Z, logS = logS, S = S)
}

##------------------------------------------------------------------------------
## Approximate log-likelihood via grid (Kitagawa discretization)
## theta.star = c(qlogis(phi), log(sigma_state), mu, log(sigma))
##   - latent: Z_t ~ N(phi*Z_{t-1}, sigma_state^2)
##   - emission: r_t | Z_t=z ~ N( (mu*e^z - 0.5 (sigma*e^z)^2) dt, (sigma*e^z)^2 dt )
##------------------------------------------------------------------------------
mllk_SSM_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  # Unpack parameters
  phi_state   <- plogis(theta.star[1])          # in (0,1)
  sigma_state <- exp(theta.star[2])             # > 0
  mu          <- theta.star[3]                  # real
  sigma       <- exp(theta.star[4])             # > 0
  
  # Grid over latent state
  b      <- seq(-bm, bm, length = m + 1)        # boundaries
  h      <- b[2] - b[1]
  bstar  <- 0.5 * (b[-1] + b[-(m + 1)])         # midpoints (m of them)
  
  # Transition matrix for latent AR(1): Z_t | Z_{t-1}=b_i ~ N(phi*b_i, sigma_state^2)
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  # Stationary init for AR(1): N(0, sigma_state^2 / (1 - phi^2))
  delta <- h * dnorm(bstar, mean = 0, sd = sigma_state / sqrt(1 - phi_state^2))
  
  # Emission density for each time given grid point z = bstar
  # r_t ~ N( (mu*e^z - 0.5 (sigma*e^z)^2) dt, (sigma*e^z)^2 dt )
  mean_vec <- (mu * exp(bstar) - 0.5 * (sigma * exp(bstar))^2) * dt
  sd_vec   <- (sigma * exp(bstar)) * sqrt(dt)
  
  # Forward algorithm (scaled)
  # t = 1
  dens1 <- dnorm(x[1], mean = mean_vec, sd = sd_vec)
  foo   <- delta * dens1
  l     <- log(sum(foo))
  phi   <- foo / sum(foo)
  
  # t = 2..T
  Tn <- length(x)
  for (t in 2:Tn) {
    denst <- dnorm(x[t], mean = mean_vec, sd = sd_vec)
    foo   <- as.numeric(phi %*% Gamma) * denst
    st    <- sum(foo)
    l   <- l + log(st)
    phi <- foo / st
  }
  
  return(-l)  # minimize
}

##------------------------------------------------------------------------------
## Set initial guess on transformed scale
theta.start <- c(
  qlogis(phi_s_true),     # logit(phi)
  log(sigma_s_true),      # log(sigma_state)
  mu_true,                # mu (untransformed)
  log(sigma_true)         # log(sigma)
)

##------------------------------------------------------------------------------
## Run multiple simulations and fit
num_simulations <- 1
m_grid          <- 200
bm_grid         <- 3

fitted_params_SSM_BS <- matrix(NA_real_, nrow = num_simulations, ncol = length(theta.start))

for (i in 1:num_simulations) {
  set.seed(i)
  sim <- simulate_SSM_BS(n, mu_true, sigma_true, phi_s_true, sigma_s_true, S0, dt)
  r   <- sim$r  # log-returns
  
  fit <- nlm(
    f = mllk_SSM_BS,
    p = theta.start,
    x = r,
    m = m_grid,
    bm = bm_grid,
    dt = dt,
    print.level = 0
  )
  fitted_params_SSM_BS[i, ] <- fit$estimate
}

## Back-transform to original scales
fitted_params_original_SSM_BS <- cbind(
  phi   = plogis(fitted_params_SSM_BS[, 1]),
  sigZ  = exp(fitted_params_SSM_BS[, 2]),
  mu    = fitted_params_SSM_BS[, 3],
  sigma = exp(fitted_params_SSM_BS[, 4])
)

print(head(fitted_params_original_SSM_BS, 10))

## A simple summary (means across simulations)
SSM_BS_means <- colMeans(fitted_params_original_SSM_BS)
print(SSM_BS_means)
