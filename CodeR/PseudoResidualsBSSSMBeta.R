# BS–SSM_beta: Pseudo-residuals

# ---- Paths, seed, data ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

load("sp500_clean.Rdata")                   # spx
load("fitted_params_BSSSM_2019_beta.RData") # models_BSSSM_2019_beta

# ---- Model & data setup ----
dt <- 1/252
x  <- as.numeric(na.omit(spx$Return))

m_grid  <- 200
bm_grid <- 3

# ---- Fit object ----
if (!is.null(models_BSSSM_2019_beta$fit_object)) {
  mod <- models_BSSSM_2019_beta$fit_object
} else {
  mod <- models_BSSSM_2019_beta[[1]]
}

# ---- Parameters ----
extract_parameters <- function(theta.star) {
  phi_state   <- plogis(theta.star[1])
  sigma_state <- exp(theta.star[2])
  mu0         <- theta.star[3]
  sigma       <- exp(theta.star[4])
  beta_mu     <- theta.star[5]
  beta_sigma  <- theta.star[6]
  list(
    mu0         = mu0,
    sigma       = sigma,
    phi_state   = phi_state,
    sigma_state = sigma_state,
    beta_mu     = beta_mu,
    beta_sigma  = beta_sigma
  )
}

# ---- Log-forward (scaled) ----
lForward <- function(y, mod, m, bm){
  params <- extract_parameters(mod$estimate)
  T      <- length(y)
  
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(
      bstar,
      mean = params$phi_state * bstar[i],
      sd   = params$sigma_state
    )
  }
  
  sd0   <- params$sigma_state / sqrt(1 - params$phi_state^2)
  delta <- h * dnorm(bstar, 0, sd0)
  
  mu_state      <- params$mu0 + params$beta_mu    * bstar
  sigma_state_v <- params$sigma * exp(params$beta_sigma * bstar)
  
  mu_grid  <- (mu_state - 0.5 * sigma_state_v^2) * dt
  sd_grid  <- sigma_state_v * sqrt(dt)
  
  mus    <- matrix(mu_grid, nrow = T, ncol = m, byrow = TRUE)
  sigmas <- matrix(sd_grid, nrow = T, ncol = m, byrow = TRUE)
  
  lalpha <- matrix(NA_real_, m, T)
  P   <- dnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  foo <- delta * P
  s1  <- sum(foo)
  lscale <- log(s1)
  foo <- foo / s1
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:T) {
    P   <- dnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    foo <- foo %*% Gamma * P
    st  <- sum(foo)
    lscale <- lscale + log(st)
    foo <- foo / st
    lalpha[, i] <- log(foo) + lscale
  }
  lalpha
}

# ---- One-step-ahead pseudo-residuals ----
PseudoResiduals <- function(y, mod, m, bm) {
  params <- extract_parameters(mod$estimate)
  n      <- length(y)
  
  b     <- seq(-bm, bm, length = m + 1)
  h     <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(
      bstar,
      mean = params$phi_state * bstar[i],
      sd   = params$sigma_state
    )
  }
  sd0   <- params$sigma_state / sqrt(1 - params$phi_state^2)
  delta <- h * dnorm(bstar, 0, sd0)
  
  mu_state      <- params$mu0 + params$beta_mu    * bstar
  sigma_state_v <- params$sigma * exp(params$beta_sigma * bstar)
  
  mu_grid  <- (mu_state - 0.5 * sigma_state_v^2) * dt
  sd_grid  <- sigma_state_v * sqrt(dt)
  
  mus    <- matrix(mu_grid, nrow = n, ncol = m, byrow = TRUE)
  sigmas <- matrix(sd_grid, nrow = n, ncol = m, byrow = TRUE)
  
  la <- t(lForward(y = y, mod = mod, m = m, bm = bm))
  Res  <- rep(NA_real_, n)
  pMat <- matrix(NA_real_, nrow = n, ncol = m)
  
  pMat[1, ] <- pnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  Res[1]    <- qnorm(as.numeric(delta %*% pMat[1, ]))
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    
    c <- max(la[i - 1, ])
    a <- exp(la[i - 1, ] - c)
    
    weighted_Gamma <- Gamma / sum(a)
    
    Res[i] <- qnorm(as.numeric(a %*% weighted_Gamma %*% pMat[i, ]))
  }
  
  list(Res = Res)
}

# ---- Compute pseudo-residuals ----
pseudo_res_SSM <- PseudoResiduals(y = x, mod = mod, m = m_grid, bm = bm_grid)

# ---- Diagnostics (2×2 panel) ----
res_SSM <- pseudo_res_SSM$Res

finite_idx     <- is.finite(res_SSM)
res_SSM_finite <- res_SSM[finite_idx]

date_all         <- spx$Date[!is.na(spx$Return)]
date_aligned_SSM <- date_all[finite_idx]

png("BSSSM_PseudoResidualDiagnostics.png", width = 2400, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.8, 0))

col_main <- "#901a1E"
col_ref  <- "steelblue"

# (1,1) Time series
if (length(res_SSM_finite) > 0) {
  if (length(date_aligned_SSM) == length(res_SSM_finite)) {
    plot(date_aligned_SSM, res_SSM_finite, type = "l",
         col = col_main, lwd = 1.2,
         xlab = "", ylab = "Pseudo-residuals")
  } else {
    plot(seq_along(res_SSM_finite), res_SSM_finite, type = "l",
         col = col_main, lwd = 1.2,
         xlab = "", ylab = "Pseudo-residuals")
  }
  abline(h = 0, lwd = 1, col = "grey60")
  abline(h = c(-3, 3), lty = 2, col = "grey75")
  title("Time series", cex.main = 1.4)
} else {
  plot(0, 0, type = "n", xlab = "", ylab = "", main = "Time series")
  text(0, 0, "No finite residuals", cex = 1.2)
}

# (1,2) Histogram + N(0,1)
if (length(res_SSM_finite) > 0) {
  hist(res_SSM_finite,
       breaks = 40,
       freq   = FALSE,
       col    = col_main,
       border = "white",
       main   = "Histogram",
       xlab   = "Pseudo-residual",
       ylab   = "Density")
  curve(dnorm(x, mean = 0, sd = 1),
        add = TRUE, col = col_ref, lwd = 2)
} else {
  plot(0, 0, type = "n", xlab = "", ylab = "", main = "Histogram")
  text(0, 0, "No finite residuals", cex = 1.2)
}

# (2,1) QQ-plot
if (length(res_SSM_finite) > 0) {
  qqnorm(res_SSM_finite,
         main = "QQ-plot",
         col  = col_main,
         xlab = "Theoretical quantiles",
         ylab = "Sample quantiles",
         cex.main = 1.4)
  qqline(res_SSM_finite, col = col_ref, lwd = 2)
} else {
  plot(0, 0, type = "n", xlab = "", ylab = "", main = "QQ-plot")
  text(0, 0, "No finite residuals", cex = 1.2)
}

# (2,2) ACF
if (length(res_SSM_finite) > 1) {
  acf_obj <- acf(res_SSM_finite, plot = FALSE)
  plot(acf_obj,
       main = "",
       xlab = "Lag",
       ylab = "Autocorrelation",
       col  = col_main,
       lwd  = 2)
  title("ACF", cex.main = 1.4)
} else {
  plot(0, 0, type = "n", xlab = "", ylab = "", main = "ACF")
  text(0, 0, "Not enough finite residuals", cex = 1.2)
}

mtext(
  expression(BS-SSM[beta] ~ "Pseudo-Residuals"),
  side  = 3,
  outer = TRUE,
  cex   = 1,
  font  = 2
)

dev.off()
