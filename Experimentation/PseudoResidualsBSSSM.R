####################################################################################################################
#----------------------------------------- Black–Scholes SSM: Pseudo-Residuals ------------------------------------
####################################################################################################################

### Libraries
suppressPackageStartupMessages({
  library(tidyverse)      # data wrangling
  library(dtplyr)
  library(ggplot2)
  library(gridExtra)
  library(fHMM)
  library(Quandl)
  library(dplyr)
  library(lubridate)
  library(stats)
  library(matrixcalc)
  library(RColorBrewer)
  library(latex2exp)
  library(matrixStats)
})

### Paths, seed, data
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# Your cleaned S&P 500 dataset (expects object `spx` with column `Return` = log returns)
load("sp500_clean.Rdata")    # -> spx

# Your fitted BS–SSM model list; contains list(fit_BS) from nlm()
# If you saved under a different filename, adjust the file name here.
load("fitted_params_BSSSM_2019.RData")   # -> models_BSSSM_2019

### Model & data setup
dt <- 1/252
x  <- as.numeric(na.omit(spx$Return))

# Grid for continuous-state approximation
m_grid  <- 200
bm_grid <- 3

# Pull your nlm() fit (fit_BS) from the saved list
mod <- models_BSSSM_2019[[1]]   # == fit_BS

#---------------------------------------------------------------------------------------------------
# Parameter extraction (matches your nlm parameterization)
# theta.star = mod$estimate with order:
#   [1] phi_state (logit scale)  -> plogis
#   [2] sigma_state (log scale)  -> exp
#   [3] mu (unconstrained)
#   [4] sigma (log scale)        -> exp
#---------------------------------------------------------------------------------------------------
extract_parameters <- function(theta.star) {
  phi_state   <- plogis(theta.star[1])
  sigma_state <- exp(theta.star[2])
  mu          <- theta.star[3]
  sigma       <- exp(theta.star[4])
  list(mu = mu, sigma = sigma, phi_state = phi_state, sigma_state = sigma_state)
}

#---------------------------------------------------------------------------------------------------
# Log-forward probabilities for continuous-state SSM (Kitagawa discretization)
# Emission: Black–Scholes log-return R_t ~ N( (μ - 0.5 σ^2) dt,  σ^2 dt )
#---------------------------------------------------------------------------------------------------
lForward <- function(y, mod, m, bm){
  params <- extract_parameters(mod$estimate)
  T <- length(y)
  
  # Discretize latent AR(1) state on [-bm, bm]
  b <- seq(-bm, bm, length = m + 1)        # boundaries
  h <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5     # midpoints
  
  # Transition probability matrix for latent AR(1) (approximate HMM)
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state)
  }
  
  # Stationary distribution of latent AR(1)
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2))
  
  # ---- BS emissions (same across latent grid because μ, σ are global here) ----
  mus    <- matrix((params$mu - 0.5 * params$sigma^2) * dt, nrow = T, ncol = m, byrow = TRUE)
  sigmas <- matrix(params$sigma * sqrt(dt),                 nrow = T, ncol = m, byrow = TRUE)
  
  # Forward recursion in log-space with scaling
  lalpha <- matrix(NA_real_, m, T)
  P   <- dnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  foo <- delta * P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:T) {
    P   <- dnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    foo <- foo %*% Gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  lalpha
}

#---------------------------------------------------------------------------------------------------
# One-step-ahead forecast pseudo-residuals
#---------------------------------------------------------------------------------------------------
PseudoResiduals <- function(y, mod, m, bm) {
  params <- extract_parameters(mod$estimate)
  
  # Discretize latent AR(1)
  b <- seq(-bm, bm, length = m + 1)
  h <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5
  
  # Transition matrix and stationary dist
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state)
  }
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2))
  
  # ---- BS emissions ----
  mus    <- matrix((params$mu - 0.5 * params$sigma^2) * dt, nrow = length(y), ncol = m, byrow = TRUE)
  sigmas <- matrix(params$sigma * sqrt(dt),                 nrow = length(y), ncol = m, byrow = TRUE)
  
  la <- t(lForward(y = y, mod = mod, m = m, bm = bm))
  n <- length(y)
  Res  <- rep(NA_real_, n)
  pMat <- matrix(NA_real_, nrow = n, ncol = m)
  
  pMat[1, ] <- pnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  Res[1]    <- qnorm(as.numeric(delta %*% pMat[1, ]))
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    c <- max(la[i - 1, ])
    a <- exp(la[i - 1, ] - c)
    weighted_Gamma <- Gamma / sum(a)  # retain your normalization trick
    Res[i] <- qnorm(as.numeric(a %*% weighted_Gamma %*% pMat[i, ]))
  }
  
  list(Res = Res)
}



################################################################################################
#---------------- Diagnostics for BS–SSM pseudo-residuals (2×2 panel) -------------------------#
################################################################################################

# Residuals
res_SSM        <- pseudo_res_SSM$Res
res_SSM_finite <- res_SSM[is.finite(res_SSM)]

# Optional: aligned dates (adjust column name if needed)
date_aligned_SSM <- spx$Date[!is.na(spx$Return)]  # change "Date" if different

png("BSSSM_PseudoResidualDiagnostics.png", width = 2400, height = 2400, res = 300)

# Slightly larger top margin so titles don't get clipped
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.8, 0))

col_main <- "#901a1E"
col_ref  <- "steelblue"

### (1,1) Time series of pseudo-residuals
if (length(date_aligned_SSM) == length(res_SSM)) {
  plot(date_aligned_SSM, res_SSM, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Pseudo-residuals")
} else {
  plot(seq_along(res_SSM), res_SSM, type = "l",
       col = col_main, lwd = 1.2,
       xlab = "", ylab = "Pseudo-residuals")
}
abline(h = 0, lwd = 1, col = "grey60")
abline(h = c(-3, 3), lty = 2, col = "grey75")
title("Time series", cex.main = 1.4)

### (1,2) Histogram (trimmed 1–99%) + N(0,1) overlay
q_low  <- quantile(res_SSM_finite, 0.01, na.rm = TRUE)
q_high <- quantile(res_SSM_finite, 0.99, na.rm = TRUE)
res_trimmed <- res_SSM_finite[res_SSM_finite >= q_low & res_SSM_finite <= q_high]

hist(res_trimmed,
     breaks = 40,
     freq   = FALSE,
     col    = col_main,
     border = "white",
     main   = "Histogram (trimmed 1–99%)",
     xlab   = "Pseudo-residual",
     ylab   = "Density")
curve(dnorm(x, mean = 0, sd = 1),
      add = TRUE, col = col_ref, lwd = 2)

### (2,1) QQ-plot vs N(0,1)
qqnorm(res_SSM_finite,
       main = "QQ-plot",
       col  = col_main,
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles",
       cex.main = 1.4)
qqline(res_SSM_finite, col = col_ref, lwd = 2)

### (2,2) ACF of pseudo-residuals
acf_obj <- acf(res_SSM_finite, plot = FALSE)

# Let plot.acf draw the panel, then force the title explicitly
plot(acf_obj,
     main = "",                 # leave empty here
     xlab = "Lag",
     ylab = "Autocorrelation",
     col  = col_main,
     lwd  = 2)
title("ACF", cex.main = 1.4)    # <- this WILL write text above the panel

mtext("BS-SSM Pseudo-Residuals",
      side  = 3,
      outer = TRUE,
      cex   = 1,
      font  = 2)   # 2 = bold
dev.off()
