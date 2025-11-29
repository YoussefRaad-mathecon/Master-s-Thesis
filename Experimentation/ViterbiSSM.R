################################################################################################
#--------------------------- Fitting + Viterbi: Black–Scholes SSM (GBM) ------------------------
################################################################################################

### 0) Libraries & setup
suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(ggplot2)
  library(gridExtra)
})

set.seed(112233)

setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("fitted_params_BSSSM_2019.RData")  # expect models_BSSSM_2019 list with $theta_bt and/or $fit_BS (nlm())

dt      <- 1/252
m_grid  <- 200
bm_grid <- 3

################################################################################################
# 1) Load data (log-returns) -------------------------------------------------------------------
################################################################################################

# Expecting 'spx' inside this Rdata with a column 'Return' (already log-returns)
load("sp500_clean.Rdata")   # provides 'spx'

spx <- spx %>% slice_head(n = 23109)
S   <- na.omit(spx$Return)  # drop initial NA from diff/log-return construction
x   <- as.numeric(S)

# Optional date if available in your object:
date_vec <- NULL
if ("Date" %in% names(spx))     date_vec <- spx$Date[!is.na(spx$Return)]
if ("DateCont" %in% names(spx)) date_vec <- spx$DateCont[!is.na(spx$Return)]

################################################################################################
# 2) BS-SSM log-likelihood (approx HMM via grid on latent Z_t) --------------------------------
################################################################################################
# Z_t ~ AR(1): Z_t = phi * Z_{t-1} + sigma_state * eps_t, eps_t ~ N(0,1)
# Emission: x_t | Z_t=z ~ N( (mu*e^z - 0.5*(sigma*e^z)^2)*dt , (sigma*e^z)*sqrt(dt) )

mllk_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  phi_state   <- plogis(theta.star[1])  # map R -> (0,1)
  sigma_state <- exp(theta.star[2])     # map R -> (0,inf)
  mu          <- theta.star[3]          # R
  sigma       <- exp(theta.star[4])     # (0,inf)
  
  # latent grid
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])
  
  # transition kernel Gamma (rows i -> columns j)
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  rs <- rowSums(Gamma)
  if (any(rs <= 0 | !is.finite(rs))) return(1e50)
  Gamma <- Gamma / rs
  
  # stationary init for AR(1)
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  sdelta <- sum(delta)
  if (!is.finite(sdelta) || sdelta <= 0) return(1e50)
  delta <- delta / sdelta
  
  # emission parameters across grid
  expZ     <- exp(bstar)
  mean_vec <- (mu * expZ - 0.5 * (sigma * expZ)^2) * dt
  sd_vec   <- (sigma * expZ) * sqrt(dt)
  
  # forward recursion (scaled)
  dens1 <- dnorm(x[1], mean = mean_vec, sd = sd_vec)
  foo   <- delta * dens1
  s1    <- sum(foo)
  if (!is.finite(s1) || s1 <= 0) return(1e50)
  logL <- log(s1)
  phi  <- foo / s1
  
  for (t in 2:length(x)) {
    denst <- dnorm(x[t], mean = mean_vec, sd = sd_vec)
    foo   <- as.numeric(phi %*% Gamma) * denst
    st    <- sum(foo)
    if (!is.finite(st) || st <= 0) return(1e50)
    logL <- logL + log(st)
    phi  <- foo / st
  }
  
  return(-logL)
}

################################################################################################
# 3) (Optional) Fit the BS-SSM with optim ------------------------------------------------------
################################################################################################

DO_FIT <- FALSE  # set TRUE if you want to refit with 'optim'

theta.start <- c(
  qlogis(0.95),  # phi_state (transformed)
  log(0.08),     # sigma_state (transformed)
  0.30,          # mu
  log(0.08)      # sigma (transformed)
)

if (DO_FIT) {
  fit_BS <- optim(
    par     = theta.start,
    fn      = mllk_BS,
    x       = x,
    m       = m_grid,
    bm      = bm_grid,
    dt      = dt,
    method  = "BFGS",
    control = list(maxit = 1000, reltol = 1e-8),
    hessian = TRUE
  )
  
  models_BSSSM_2019 <- list(fit_BS = fit_BS)
  save(models_BSSSM_2019, file = "fitted_params_BSSSM_2019.RData")
}

################################################################################################
# 4) Parameter extraction (robust to nlm() or pre-made theta_bt) -------------------------------
################################################################################################
# Preferred: use already back-transformed theta_bt if present.
# Otherwise: back-transform nlm()$estimate.
# Otherwise: fall back to theta.start.

if (!exists("theta_BS_hat")) {
  
  # If the models list not present (e.g. script run fresh), try loading again
  if (!exists("models_BSSSM_2019") && file.exists("fitted_params_BSSSM_2019.RData")) {
    load("fitted_params_BSSSM_2019.RData")
  }
  
  if (exists("models_BSSSM_2019") && !is.null(models_BSSSM_2019$theta_bt)) {
    # Already back-transformed (phi, sigZ, mu, sigma)
    theta_BS_hat <- as.numeric(models_BSSSM_2019$theta_bt)
    names(theta_BS_hat) <- c("phi", "sigZ", "mu", "sigma")
    
  } else if (exists("models_BSSSM_2019") && !is.null(models_BSSSM_2019$fit_BS$estimate)) {
    # nlm() object: back-transform $estimate
    theta_hat_tr <- as.numeric(models_BSSSM_2019$fit_BS$estimate)
    stopifnot(length(theta_hat_tr) >= 4)
    theta_BS_hat <- c(
      phi   = plogis(theta_hat_tr[1]),
      sigZ  = exp(theta_hat_tr[2]),
      mu    = theta_hat_tr[3],
      sigma = exp(theta_hat_tr[4])
    )
    
  } else if (exists("models_BSSSM_2019") && !is.null(models_BSSSM_2019$fit_BS$par)) {
    # in case you later save an optim() fit with $par
    theta_hat_tr <- as.numeric(models_BSSSM_2019$fit_BS$par)
    stopifnot(length(theta_hat_tr) >= 4)
    theta_BS_hat <- c(
      phi   = plogis(theta_hat_tr[1]),
      sigZ  = exp(theta_hat_tr[2]),
      mu    = theta_hat_tr[3],
      sigma = exp(theta_hat_tr[4])
    )
    
  } else {
    message("Note: Saved fit not found; using theta.start as fallback.")
    theta_BS_hat <- c(
      phi   = plogis(theta.start[1]),
      sigZ  = exp(theta.start[2]),
      mu    = theta.start[3],
      sigma = exp(theta.start[4])
    )
  }
}

print(theta_BS_hat)

################################################################################################
# 5) Viterbi algorithm for BS-SSM --------------------------------------------------------------
################################################################################################

viterbi_bs <- function(x, mu, sigma, phi_state, sigma_state, m = 200, bm = 3, dt = 1/252) {
  n <- length(x)
  
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])
  
  # Transition Gamma
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  rs <- rowSums(Gamma)
  if (any(rs <= 0 | !is.finite(rs))) stop("Invalid Gamma rows in Viterbi.")
  Gamma <- Gamma / rs
  
  # Stationary init
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, 0, sd0)
  delta <- delta / sum(delta)
  
  # Emission grid
  expZ     <- exp(bstar)
  mean_vec <- (mu * expZ - 0.5 * (sigma * expZ)^2) * dt
  sd_vec   <- (sigma * expZ) * sqrt(dt)
  
  log_emiss <- function(x_t) {
    dnorm(x_t, mean = mean_vec, sd = sd_vec, log = TRUE)
  }
  
  # Log-domain Viterbi
  logGamma <- log(Gamma)
  logdelta <- log(delta)
  
  nu  <- matrix(-Inf, n, m)
  psi <- matrix(0L, n, m)
  
  nu[1, ] <- logdelta + log_emiss(x[1])
  
  for (t in 2:n) {
    le <- log_emiss(x[t])
    for (j in 1:m) {
      v <- nu[t - 1, ] + logGamma[, j]
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
  
  # Return latent Z_t midpoints
  bstar[iv]
}

################################################################################################
# 6) Run Viterbi + derive state-dependent paths ------------------------------------------------
################################################################################################

Zhat_BS <- viterbi_bs(
  x          = x,
  mu         = theta_BS_hat["mu"],
  sigma      = theta_BS_hat["sigma"],
  phi_state  = theta_BS_hat["phi"],
  sigma_state = theta_BS_hat["sigZ"],
  m          = m_grid,
  bm         = bm_grid,
  dt         = dt
)

sigma_t <- as.numeric(theta_BS_hat["sigma"] * exp(Zhat_BS))
mu_t    <- as.numeric(theta_BS_hat["mu"]    * exp(Zhat_BS))

cond_mean_t <- (mu_t - 0.5 * sigma_t^2) * dt
cond_sd_t   <- sigma_t * sqrt(dt)

bs_paths_df <- tibble(
  t           = seq_along(x),
  date        = if (!is.null(date_vec)) date_vec else seq_along(x),
  x           = x,
  Zhat        = Zhat_BS,
  sigma_t     = sigma_t,
  mu_t        = mu_t,
  cond_mean_t = cond_mean_t,
  cond_sd_t   = cond_sd_t
)

################################################################################################
# 7) Plots -------------------------------------------------------------------------------------
################################################################################################

# (A) Scatter of returns colored by sigma_t (legend to the RIGHT)
p_scatter_sigma <- ggplot(bs_paths_df, aes(x = t, y = x, color = sigma_t)) +
  geom_point(alpha = 0.6, size = 1) +
  labs(
    title = expression(paste("Returns colored by ", sigma[t])),
    x     = "t",
    y     = "Log-returns",
    color = expression(sigma[t])
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Scatter of returns colored by mu_t (legend to the RIGHT)
p_scatter_mu <- ggplot(bs_paths_df, aes(x = t, y = x, color = mu_t)) +
  geom_point(alpha = 0.6, size = 1) +
  labs(
    title = expression(paste("Returns colored by ", mu[t])),
    x     = "t",
    y     = "Log-returns",
    color = expression(mu[t])
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Print if running interactively
p_scatter_mu
p_scatter_sigma
# p_scatter  # (uncomment if you define this elsewhere)

#--------------------------------- FRONT PAGE ILLUSTRATION -------------------------------------#
################################################################################################
# 7) Front-page illustration plot (boxed, no ticks, tall color bar) ----------------------------
################################################################################################

red_palette <- c("#F7D9DC", "#E4969B", "#C35359", "#901A1E")

p_front_sigma <- ggplot(bs_paths_df, aes(x = t, y = x, color = sigma_t)) +
  geom_point(alpha = 0.8, size = 0.7) +
  scale_color_gradientn(colours = red_palette) +
  # No ticks, but keep axes & sec axes to get a full box
  scale_x_continuous(
    breaks   = NULL,
    sec.axis = dup_axis(labels = NULL, breaks = NULL)
  ) +
  scale_y_continuous(
    breaks   = NULL,
    sec.axis = dup_axis(labels = NULL, breaks = NULL)
  ) +
  guides(
    colour = guide_colourbar(
      barheight = grid::unit(0.9, "npc"),  # make the bar long (near full axis height)
      barwidth  = grid::unit(0.03, "npc")  # thin bar
    )
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    # No axis titles or text
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    # Box around plot, but no ticks
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks      = element_blank(),
    # Legend: keep bar only, no text
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_blank(),
    plot.title      = element_blank(),
    plot.subtitle   = element_blank(),
    plot.caption    = element_blank(),
    plot.margin     = margin(5, 5, 5, 5)
  )

p_front_sigma
