####################################################################################################################
#----------------------------------------- Black–Scholes N-state HMM ----------------------------------------------
####################################################################################################################

### Minimal lib
suppressPackageStartupMessages({
  library(tidyverse)
})

### Paths + data
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("sp500_clean.Rdata")                # -> spx with column `Return` (log returns)
load("fitted_params_BSHMM_2019.RData")   # -> models_BSHMM_2019 (your list of nlm fits)

dt <- 1/252
x  <- as.numeric(na.omit(spx$Return))

# --------------------------------------------------------------------------------------------------
# Parameter extraction for your three families:
#   pattern = "mu"        : mu_j state-dependent, sigma shared
#            = "sigma"    : mu shared, sigma_j state-dependent
#            = "mu_sigma" : mu_j and sigma_j both state-dependent
# The transition matrix parameterization matches your likelihood code exactly.
# --------------------------------------------------------------------------------------------------
extract_parameters_HMM <- function(theta.star, N, pattern = c("mu","sigma","mu_sigma")) {
  pattern <- match.arg(pattern)
  
  if (pattern == "mu") {
    mu_vec    <- theta.star[1:N]
    sigma_vec <- rep(exp(theta.star[N + 1]), N)
    ofs_start <- N + 2
    
  } else if (pattern == "sigma") {
    mu_vec    <- rep(theta.star[1], N)
    sigma_vec <- exp(theta.star[2:(N + 1)])
    ofs_start <- N + 2
    
  } else { # "mu_sigma"
    mu_vec    <- theta.star[1:N]
    sigma_vec <- exp(theta.star[(N + 1):(2 * N)])
    ofs_start <- 2 * N + 1
  }
  
  # Transition matrix (keep identical to your mllk construction)
  Gamma <- diag(N)
  if (ofs_start <= length(theta.star)) {
    off_diag_elements <- exp(theta.star[ofs_start:length(theta.star)])
    Gamma[!Gamma] <- off_diag_elements
  }
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution (match your likelihood: eigen of t(Gamma))
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  list(mu = mu_vec, sigma = sigma_vec, Gamma = Gamma, delta = delta)
}

# --------------------------------------------------------------------------------------------------
# Log-forward (scaled) for HMM with BS emissions
# --------------------------------------------------------------------------------------------------
lForward_HMM <- function(y, mod, N, pattern = c("mu","sigma","mu_sigma"), dt = 1/252) {
  pattern <- match.arg(pattern)
  params  <- extract_parameters_HMM(mod$estimate, N, pattern)
  Tn <- length(y)
  
  # Emission parameters per-state
  mu_t     <- (params$mu - 0.5 * (params$sigma^2)) * dt     # length N
  sigma_t  <- params$sigma * sqrt(dt)                       # length N
  mus_mat    <- matrix(mu_t,    nrow = Tn, ncol = N, byrow = TRUE)
  sigmas_mat <- matrix(sigma_t, nrow = Tn, ncol = N, byrow = TRUE)
  
  lalpha <- matrix(NA_real_, N, Tn)
  
  P1  <- dnorm(y[1], mean = mus_mat[1, ], sd = sigmas_mat[1, ])
  foo <- params$delta * P1
  s1  <- sum(foo)
  lscale <- log(s1)
  foo <- foo / s1
  lalpha[, 1] <- lscale + log(foo)
  
  for (t in 2:Tn) {
    Pt  <- dnorm(y[t], mean = mus_mat[t, ], sd = sigmas_mat[t, ])
    foo <- foo %*% params$Gamma * Pt
    st  <- sum(foo)
    lscale <- lscale + log(st)
    foo    <- foo / st
    lalpha[, t] <- log(foo) + lscale
  }
  lalpha
}

# --------------------------------------------------------------------------------------------------
# One-step-ahead forecast pseudo-residuals (keeps your weighting trick)
# --------------------------------------------------------------------------------------------------
PseudoResiduals_HMM <- function(y, mod, N, pattern = c("mu","sigma","mu_sigma"), dt = 1/252) {
  pattern <- match.arg(pattern)
  params  <- extract_parameters_HMM(mod$estimate, N, pattern)
  n <- length(y)
  
  mu_t     <- (params$mu - 0.5 * (params$sigma^2)) * dt
  sigma_t  <- params$sigma * sqrt(dt)
  mus_mat    <- matrix(mu_t,    nrow = n, ncol = N, byrow = TRUE)
  sigmas_mat <- matrix(sigma_t, nrow = n, ncol = N, byrow = TRUE)
  
  la   <- t(lForward_HMM(y = y, mod = mod, N = N, pattern = pattern, dt = dt))
  Res  <- rep(NA_real_, n)
  pMat <- matrix(NA_real_, nrow = n, ncol = N)
  
  pMat[1, ] <- pnorm(y[1], mean = mus_mat[1, ], sd = sigmas_mat[1, ])
  Res[1]    <- qnorm(as.numeric(params$delta %*% pMat[1, ]))
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus_mat[i, ], sd = sigmas_mat[i, ])
    c  <- max(la[i - 1, ])
    a  <- exp(la[i - 1, ] - c)
    weighted_Gamma <- params$Gamma / sum(a)  # your original normalization trick
    Res[i] <- qnorm(as.numeric(a %*% weighted_Gamma %*% pMat[i, ]))
  }
  
  list(Res = Res)
}

# --------------------------------------------------------------------------------------------------
# Examples: compute pseudo-residuals for each saved model (2–5 states)
# --------------------------------------------------------------------------------------------------

# μ-only models
pr_2_mu <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod2_mu, N = 2, pattern = "mu", dt = dt)
pr_3_mu <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod3_mu, N = 3, pattern = "mu", dt = dt) # 2 NAN
pr_4_mu <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod4_mu, N = 4, pattern = "mu", dt = dt)
pr_5_mu <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod5_mu, N = 5, pattern = "mu", dt = dt)

# σ-only models
pr_2_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod2_sigma, N = 2, pattern = "sigma", dt = dt)
pr_3_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod3_sigma, N = 3, pattern = "sigma", dt = dt)
pr_4_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod4_sigma, N = 4, pattern = "sigma", dt = dt)
pr_5_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod5_sigma, N = 5, pattern = "sigma", dt = dt)

# μ & σ models
pr_2_mu_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod2_mu_sigma, N = 2, pattern = "mu_sigma", dt = dt)
pr_3_mu_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod3_mu_sigma, N = 3, pattern = "mu_sigma", dt = dt)
pr_4_mu_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod4_mu_sigma, N = 4, pattern = "mu_sigma", dt = dt)
pr_5_mu_sigma <- PseudoResiduals_HMM(x, models_BSHMM_2019$mod5_mu_sigma, N = 5, pattern = "mu_sigma", dt = dt)


#----------------------------------------- SMALL PLOTS (μ, σ, μ&σ) -----------------------------------------------
par(oma = c(3,3,0,1), mar = c(3,3,2,3), mfrow = c(3,4))

## Row 1: mu-only
clean_res <- function(z) z[is.finite(z)]

## Row 1: mu-only
qqnorm(clean_res(pr_2_mu$Res), main = "2-State HMM",
       col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7)
qqline(clean_res(pr_2_mu$Res), col = "steelblue", lwd = 3)

qqnorm(clean_res(pr_3_mu$Res), main = "3-State HMM",
       col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7)
qqline(clean_res(pr_3_mu$Res), col = "steelblue", lwd = 3)

qqnorm(clean_res(pr_4_mu$Res), main = "4-State HMM",
       col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7)
qqline(clean_res(pr_4_mu$Res), col = "steelblue", lwd = 3)

qqnorm(clean_res(pr_5_mu$Res), main = "5-State HMM",
       col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7)
qqline(clean_res(pr_5_mu$Res), col = "steelblue", lwd = 3)

## Row 2: sigma-only
qqnorm(pr_2_sigma$Res, main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_2_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_3_sigma$Res, main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_3_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_4_sigma$Res, main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_4_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_5_sigma$Res, main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_5_sigma$Res, col = "steelblue", lwd = 3)

## Row 3: mu & sigma
qqnorm(pr_2_mu_sigma$Res, main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_2_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_3_mu_sigma$Res, main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_3_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_4_mu_sigma$Res, main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_4_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_5_mu_sigma$Res, main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_5_mu_sigma$Res, col = "steelblue", lwd = 3)

mtext(text = "Theoretical Quantiles", side = 1, line = 1, outer = TRUE, cex = 2.7)
mtext(text = "Sample Quantiles",     side = 2, line = 0, outer = TRUE, cex = 2.7, adj = 0.5)
mtext(text = expression(mu),         side = 4, line = -1, outer = TRUE, cex = 2, adj = 0.85)
mtext(text = expression(sigma),      side = 4, line = -1, outer = TRUE, cex = 2, adj = 0.50)
mtext(text = expression(mu * ", " * sigma), side = 4, line = -1, outer = TRUE, cex = 2, adj = 0.15)

#----------------------------------------- BIG PLOT (μ, σ, μ&σ) --------------------------------------------------
par(mfrow = c(3,4), mar = c(2,2,2,1), oma = c(3,3,0,2))

## mu-only
qqnorm(pr_2_mu$Res, main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_2_mu$Res, col = "steelblue", lwd = 3)
qqnorm(pr_3_mu$Res, main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_3_mu$Res, col = "steelblue", lwd = 3)
qqnorm(pr_4_mu$Res, main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_4_mu$Res, col = "steelblue", lwd = 3)
qqnorm(pr_5_mu$Res, main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_5_mu$Res, col = "steelblue", lwd = 3)

## sigma-only
qqnorm(pr_2_sigma$Res, main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_2_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_3_sigma$Res, main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_3_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_4_sigma$Res, main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_4_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_5_sigma$Res, main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_5_sigma$Res, col = "steelblue", lwd = 3)

## mu & sigma
qqnorm(pr_2_mu_sigma$Res, main = "2-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_2_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_3_mu_sigma$Res, main = "3-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_3_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_4_mu_sigma$Res, main = "4-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_4_mu_sigma$Res, col = "steelblue", lwd = 3)
qqnorm(pr_5_mu_sigma$Res, main = "5-State HMM", col = "#901a1E", xlab = "", ylab = "", cex.main = 2.7); qqline(pr_5_mu_sigma$Res, col = "steelblue", lwd = 3)

mtext(text = "Theoretical Quantiles", side = 1, line = 2, outer = TRUE, cex = 2.7)
mtext(text = "Sample Quantiles",      side = 2, line = 0.3, outer = TRUE, cex = 2.7, adj = 0.5)
mtext(text = expression(mu),                 side = 4, line = 0.3, outer = TRUE, cex = 2, adj = 0.87)
mtext(text = expression(sigma),              side = 4, line = -0.1, outer = TRUE, cex = 2, adj = 0.63)
mtext(text = expression(mu * ", " * sigma),  side = 4, line = -0.1, outer = TRUE, cex = 2, adj = 0.37)



################################################################################################
#-------------------------- Collect BS–HMM pseudo-residuals (4×3 grid) ------------------------#
################################################################################################

clean_res <- function(z) z[is.finite(z)]

# Pack residuals by pattern and number of states
residuals_HMM <- list(
  mu = list(
    `2` = pr_2_mu$Res,
    `3` = pr_3_mu$Res,
    `4` = pr_4_mu$Res,
    `5` = pr_5_mu$Res
  ),
  sigma = list(
    `2` = pr_2_sigma$Res,
    `3` = pr_3_sigma$Res,
    `4` = pr_4_sigma$Res,
    `5` = pr_5_sigma$Res
  ),
  mu_sigma = list(
    `2` = pr_2_mu_sigma$Res,
    `3` = pr_3_mu_sigma$Res,
    `4` = pr_4_mu_sigma$Res,
    `5` = pr_5_mu_sigma$Res
  )
)

patterns        <- c("mu", "sigma", "mu_sigma")
pattern_labels  <- list(
  expression(mu),
  expression(sigma),
  expression(mu * ", " * sigma)
)
Ns              <- 2:5
state_labels    <- paste0(Ns, "-state HMM")

# Aligned dates for time series (adjust column name if needed)
date_aligned <- spx$Date[!is.na(spx$Return)]  # change "Date" if your column is named differently

col_main <- "#901a1E"
col_ref  <- "steelblue"

# Helpers for outer row/column label positions
row_at <- rev(((1:length(Ns))     - 0.5) / length(Ns))   # side = 2
col_at <-      ((1:length(patterns)) - 0.5) / length(patterns)  # side = 3


################################################################################################
#-------------------------- 4×3 BS–HMM pseudo-residual diagnostics ---------------------------#
#   Rows   : 2-, 3-, 4-, 5-state HMM
#   Columns: mu state-dep, sigma state-dep, mu & sigma state-dep
#   Figures: time series, histograms, QQ-plots, ACF
################################################################################################

clean_res <- function(z) z[is.finite(z)]

# Pack residuals by pattern and number of states
residuals_HMM <- list(
  mu = list(
    `2` = pr_2_mu$Res,
    `3` = pr_3_mu$Res,
    `4` = pr_4_mu$Res,
    `5` = pr_5_mu$Res
  ),
  sigma = list(
    `2` = pr_2_sigma$Res,
    `3` = pr_3_sigma$Res,
    `4` = pr_4_sigma$Res,
    `5` = pr_5_sigma$Res
  ),
  mu_sigma = list(
    `2` = pr_2_mu_sigma$Res,
    `3` = pr_3_mu_sigma$Res,
    `4` = pr_4_mu_sigma$Res,
    `5` = pr_5_mu_sigma$Res
  )
)

patterns       <- c("mu", "sigma", "mu_sigma")
pattern_labels <- list(
  expression(mu),
  expression(sigma),
  expression(mu * ", " * sigma)
)
Ns           <- 2:5
state_labels <- paste0(Ns, "-state HMM")

# Positions for row (state) and column (pattern) labels
row_at <- rev(((1:length(Ns))       - 0.5) / length(Ns))       # side = 2
col_at <-      ((1:length(patterns)) - 0.5) / length(patterns) # side = 3

# Aligned dates for time series (adjust column name if needed)
date_aligned <- spx$Date[!is.na(spx$Return)]  # change "Date" if your column is named differently

col_main <- "#901a1E"
col_ref  <- "steelblue"


################################################################################################
# 1) Time series figure (4×3) -----------------------------------------------------------------#
################################################################################################

png("BSHMM_TimeSeries.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(4, 3),
    mar  = c(2, 2, 1.5, 1),
    oma  = c(4, 4, 4, 3),   # more top space
    font.main = 2)

for (i in seq_along(Ns)) {
  N <- Ns[i]
  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    res <- residuals_HMM[[pat]][[as.character(N)]]
    
    if (is.null(res) || sum(is.finite(res)) < 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      next
    }
    
    idx       <- is.finite(res)
    res_clean <- res[idx]
    
    if (length(date_aligned) == length(res)) {
      dates_clean <- date_aligned[idx]
      plot(dates_clean, res_clean, type = "l",
           col = col_main, lwd = 1,
           xlab = "", ylab = "")
    } else {
      plot(which(idx), res_clean, type = "l",
           col = col_main, lwd = 1,
           xlab = "", ylab = "")
    }
    
    abline(h = 0, lwd = 1, col = "grey60")
    abline(h = c(-3, 3), lty = 2, col = "grey80")
    title(state_labels[i], cex.main = 1.1)
  }
}

mtext("Time",             side = 1, line = 2.2, outer = TRUE, cex = 1.3)
mtext("Pseudo-residuals", side = 2, line = 2.2, outer = TRUE, cex = 1.3)

for (i in seq_along(Ns)) {
  mtext(state_labels[i], side = 2, at = row_at[i],
        line = 4.7, outer = TRUE, cex = 1.3)
}
for (j in seq_along(patterns)) {
  mtext(pattern_labels[[j]], side = 3, at = col_at[j],
        line = 0.6, outer = TRUE, cex = 1.5)
}

mtext("BS–HMM Pseudo-Residuals: Time-Series",
      side = 3, line = 2.2, outer = TRUE, cex = 1, font = 2)

dev.off()


################################################################################################
# 2) Histogram figure (4×3) -------------------------------------------------------------------#
################################################################################################

png("BSHMM_Histograms.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(4, 3),
    mar  = c(2, 2, 1.5, 1),
    oma  = c(4, 4, 4, 3),
    font.main = 2)

for (i in seq_along(Ns)) {
  N <- Ns[i]
  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    res <- residuals_HMM[[pat]][[as.character(N)]]
    
    res_clean <- clean_res(res)
    if (length(res_clean) < 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      next
    }
    
    # Trim to 1–99% to avoid extreme tails dominating
    q_low  <- quantile(res_clean, 0.01, na.rm = TRUE)
    q_high <- quantile(res_clean, 0.99, na.rm = TRUE)
    res_trim <- res_clean
    
    hist(res_trim,
         breaks = 30,
         freq   = FALSE,
         col    = col_main,
         border = "white",
         main   = "",
         xlab   = "",
         ylab   = "")
    
    curve(dnorm(x, 0, 1), add = TRUE, col = col_ref, lwd = 2)
    title(state_labels[i], cex.main = 1.1)
  }
}

mtext("Pseudo-residual", side = 1, line = 2.2, outer = TRUE, cex = 1.3)
mtext("Density",         side = 2, line = 2.2, outer = TRUE, cex = 1.3)

for (i in seq_along(Ns)) {
  mtext(state_labels[i], side = 2, at = row_at[i],
        line = 4.7, outer = TRUE, cex = 1.3)
}
for (j in seq_along(patterns)) {
  mtext(pattern_labels[[j]], side = 3, at = col_at[j],
        line = 0.6, outer = TRUE, cex = 1.5)
}

mtext("BS–HMM Pseudo-Residuals: Histograms",
      side = 3, line = 2.2, outer = TRUE, cex = 1, font = 2)

dev.off()


################################################################################################
# 3) QQ-plot figure (4×3) ---------------------------------------------------------------------#
################################################################################################

png("BSHMM_QQplots.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(4, 3),
    mar  = c(2, 2, 1.5, 1),
    oma  = c(4, 4, 4, 3),
    font.main = 2)

for (i in seq_along(Ns)) {
  N <- Ns[i]
  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    res <- residuals_HMM[[pat]][[as.character(N)]]
    
    res_clean <- clean_res(res)
    if (length(res_clean) < 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      next
    }
    
    qqnorm(res_clean,
           main = "",
           col  = col_main,
           xlab = "",
           ylab = "")
    qqline(res_clean, col = col_ref, lwd = 2)
    title(state_labels[i], cex.main = 1.1)
  }
}

mtext("Theoretical quantiles", side = 1, line = 2.2, outer = TRUE, cex = 1.3)
mtext("Sample quantiles",      side = 2, line = 2.2, outer = TRUE, cex = 1.3)

for (i in seq_along(Ns)) {
  mtext(state_labels[i], side = 2, at = row_at[i],
        line = 4.7, outer = TRUE, cex = 1.3)
}
for (j in seq_along(patterns)) {
  mtext(pattern_labels[[j]], side = 3, at = col_at[j],
        line = 0.6, outer = TRUE, cex = 1.5)
}

mtext("BS–HMM Pseudo-Residuals: QQ-plots",
      side = 3, line = 2.2, outer = TRUE, cex = 1, font = 2)

dev.off()


################################################################################################
# 4) ACF figure (4×3) -------------------------------------------------------------------------#
################################################################################################

png("BSHMM_ACF.png", width = 2400, height = 3200, res = 300)
par(mfrow = c(4, 3),
    mar  = c(2, 2, 1.5, 1),
    oma  = c(4, 4, 4, 3),
    font.main = 2)

for (i in seq_along(Ns)) {
  N <- Ns[i]
  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    res <- residuals_HMM[[pat]][[as.character(N)]]
    
    res_clean <- clean_res(res)
    if (length(res_clean) < 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      next
    }
    
    acf_obj <- acf(res_clean, plot = FALSE)
    plot(acf_obj,
         main = "",
         xlab = "",
         ylab = "",
         col  = col_main,
         lwd  = 2)
    title(state_labels[i], cex.main = 1.1)
  }
}

mtext("Lag",             side = 1, line = 2.2, outer = TRUE, cex = 1.3)
mtext("Autocorrelation", side = 2, line = 2.2, outer = TRUE, cex = 1.3)

for (i in seq_along(Ns)) {
  mtext(state_labels[i], side = 2, at = row_at[i],
        line = 4.7, outer = TRUE, cex = 1.3)
}
for (j in seq_along(patterns)) {
  mtext(pattern_labels[[j]], side = 3, at = col_at[j],
        line = 0.6, outer = TRUE, cex = 1.5)
}

mtext("BS–HMM Pseudo-Residuals: ACF",
      side = 3, line = 2.2, outer = TRUE, cex = 1, font = 2)

dev.off()

