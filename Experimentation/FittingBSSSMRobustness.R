################################################################################################
#--------------------------- Fitting SSM: Black–Scholes (GBM) ---------------------------------
################################################################################################

### Libraries
library(tidyverse)     # Data manipulation
library(dtplyr)
library(ggplot2)
library(gridExtra)
library(fHMM)          # (not used below but you had it loaded)
library(Quandl)
library(dplyr)
library(lubridate)
library(stats)
library(matrixcalc)
library(RColorBrewer)
library(latex2exp)
library(matrixStats)

### Working dir, seed, data
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)
load("sp500_clean.Rdata")

spx   <- spx %>% slice_head(n = 23109)
S     <- na.omit(spx$Return)     # already log returns
x_raw <- S
x     <- as.numeric(x_raw)       # numeric vector

### Parameters for the BS-SSM
dt      <- 1/252
# ---------------------------------------------------------------------------------------------#
# Likelihood for BS-SSM (state-space with latent Z_t AR(1)-like and Normal emissions on returns)
# ---------------------------------------------------------------------------------------------#
mllk_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  # Transformed parameters
  phi_state   <- plogis(theta.star[1])      # phi in (0,1)
  sigma_state <- exp(theta.star[2])         # >0
  mu          <- theta.star[3]              # real
  sigma       <- exp(theta.star[4])         # >0
  
  # Latent state grid for Z_t
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- 0.5 * (b[-1] + b[-(m + 1)])      # midpoints
  
  # Transition kernel for Z_t (discretized Gaussian AR(1)-style)
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  # Stationary initial distribution for Z_t
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  
  # Emission (returns conditional on Z_t)
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
  
  return(-logL)  # nlm minimizes
}

# -------------------------------- Initial guess --------------------------------------------------------#
theta.start <- c(
  qlogis(0.95),  # phi
  log(0.08),     # sigma_state
  0.3,           # mu
  log(0.08)      # sigma
)

# -------------------------------------------- Base fit -------------------------------------------------#

################################################################################################
#--------------------------- Robustness grid: 20 (m, bm) combinations --------------------------
################################################################################################

m_values  <- c(20, 40, 70, 100, 200)
bm_values <- c(0.5, 1, 2, 3, 4)
combo_grid <- expand.grid(m = m_values, bm = bm_values, KEEP.OUT.ATTRS = FALSE)

# Containers
models_BSSSM_grid <- vector("list", nrow(combo_grid))
names(models_BSSSM_grid) <- sprintf("m%g_bm%g", combo_grid$m, combo_grid$bm)

grid_summary <- tibble::tibble(
  name       = names(models_BSSSM_grid),
  m          = combo_grid$m,
  bm         = combo_grid$bm,
  status     = NA_character_,
  conv_code  = NA_integer_,      # nlm() termination code
  negloglik  = NA_real_,         # nlm()$minimum
  phi        = NA_real_,
  sigZ       = NA_real_,
  mu         = NA_real_,
  sigma      = NA_real_,
  elapsed_s  = NA_real_,
  message    = NA_character_
)

# Back-transform helper
bt <- function(est) c(
  phi   = plogis(est[1]),
  sigZ  = exp(est[2]),
  mu    = est[3],
  sigma = exp(est[4])
)

# Loop over all combos
for (i in seq_len(nrow(combo_grid))) {
  mi  <- combo_grid$m[i]
  bmi <- combo_grid$bm[i]
  t0  <- Sys.time()
  cat(sprintf("Fitting combo %2d/%2d: m = %-3d, bm = %-3.1f ... ",
              i, nrow(combo_grid), mi, bmi))
  
  fit_i <- try(
    nlm(
      f = mllk_BS,
      p = theta.start,
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
    grid_summary$status[i]    <- "fail"
    grid_summary$message[i]   <- msg
    grid_summary$elapsed_s[i] <- elapsed
    
    models_BSSSM_grid[[i]] <- list(
      status     = "fail",
      m          = mi,
      bm         = bmi,
      error      = msg,
      elapsed_s  = elapsed
    )
    cat("FAIL\n")
    next
  }
  
  # Successful nlm() call (regardless of convergence code)
  theta_hat_tr_i <- fit_i$estimate
  theta_hat_i    <- bt(theta_hat_tr_i)
  
  grid_summary$status[i]     <- "ok"
  grid_summary$conv_code[i]  <- fit_i$code
  grid_summary$negloglik[i]  <- fit_i$minimum
  grid_summary$phi[i]        <- theta_hat_i["phi"]
  grid_summary$sigZ[i]       <- theta_hat_i["sigZ"]
  grid_summary$mu[i]         <- theta_hat_i["mu"]
  grid_summary$sigma[i]      <- theta_hat_i["sigma"]
  grid_summary$elapsed_s[i]  <- elapsed
  
  models_BSSSM_grid[[i]] <- list(
    status        = "ok",
    m             = mi,
    bm            = bmi,
    fit           = fit_i,               # raw nlm object
    theta_tr      = theta_hat_tr_i,      # transformed params
    theta_bt      = theta_hat_i,         # back-transformed params
    negloglik     = fit_i$minimum,
    conv_code     = fit_i$code,
    hessian       = fit_i$hessian,
    elapsed_s     = elapsed
  )
  
  cat(sprintf("OK (code=%d, -logL=%.6f, %.2fs)\n", fit_i$code, fit_i$minimum, elapsed))
}

# Save everything from robustness run
models_BSSSM_robustness <- list(
  grid   = models_BSSSM_grid,
  summary= grid_summary
)
print(dplyr::arrange(grid_summary, bm, m))



# Optional quick inspection
grid_summary_sorted <- dplyr::arrange(grid_summary, bm, m)
save(
  models_BSSSM_robustness,   # full objects
  grid_summary,              # raw tibble
  grid_summary_sorted,       # pre-sorted tibble
  file = "fitted_params_BSSSM_grid_2019.RData"
)

# Save results
saveRDS(grid_summary,         "BSSSM_grid_summary_2019.rds")
saveRDS(grid_summary_sorted,  "BSSSM_grid_summary_sorted_2019.rds")


# Load the data and print results
# Load the data and print results 
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation") 
load("fitted_params_BSSSM_grid_2019.RData"); print(grid_summary_sorted, n = Inf) 
print(readRDS("BSSSM_grid_summary_sorted_2019.rds"), n = Inf)


## Sort by negative log-likelihood (smallest = best) and show full digits
options(digits = 22, scipen = 999, pillar.sigfig = 17, tibble.width = Inf)

fmt <- function(x, k = 16) formatC(x, format = "f", digits = k, big.mark = ",", decimal.mark = ".")
grid_summary_sorted %>%
  arrange(negloglik) %>%
  mutate(negloglik_full = fmt(negloglik)) %>%
  select(name, m, bm, conv_code, negloglik_full) %>%
  print(n = Inf, width = Inf)




