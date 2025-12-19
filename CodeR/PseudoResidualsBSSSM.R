# BS–SSM: Pseudo-residuals

# ---- Paths, seed, data ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

load("sp500_clean.Rdata")              # spx
load("fitted_params_BSSSM_2019.RData") # models_BSSSM_2019

# ---- Model & data setup ----
dt <- 1/252
x  <- as.numeric(na.omit(spx$Return))

m_grid  <- 70
bm_grid <- 3

mod <- models_BSSSM_2019[[1]]

# ---- Parameters ----
extract_parameters <- function(theta.star) {
  phi_state   <- plogis(theta.star[1])
  sigma_state <- exp(theta.star[2])
  mu          <- theta.star[3]
  sigma       <- exp(theta.star[4])
  list(mu = mu, sigma = sigma, phi_state = phi_state, sigma_state = sigma_state)
}

# ---- Log-forward (scaled) ----
lForward <- function(y, mod, m, bm){
  params <- extract_parameters(mod$estimate)
  T <- length(y)
  
  b <- seq(-bm, bm, length = m + 1)
  h <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state)
  }
  
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2))
  
  mus    <- matrix((params$mu - 0.5 * params$sigma^2) * dt, nrow = T, ncol = m, byrow = TRUE)
  sigmas <- matrix(params$sigma * sqrt(dt),                 nrow = T, ncol = m, byrow = TRUE)
  
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

# ---- One-step-ahead pseudo-residuals ----
PseudoResiduals <- function(y, mod, m, bm) {
  params <- extract_parameters(mod$estimate)
  
  b <- seq(-bm, bm, length = m + 1)
  h <- b[2] - b[1]
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state)
  }
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2))
  
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
    weighted_Gamma <- Gamma / sum(a)
    Res[i] <- qnorm(as.numeric(a %*% weighted_Gamma %*% pMat[i, ]))
  }
  
  list(Res = Res)
}

# ---- Diagnostics (2×2 panel) ----
pseudo_res_SSM <- PseudoResiduals(
  y   = x,
  mod = mod,
  m   = m_grid,
  bm  = bm_grid
)

res_SSM        <- pseudo_res_SSM$Res
res_SSM_finite <- res_SSM[is.finite(res_SSM)]

date_aligned_SSM <- spx$Date[!is.na(spx$Return)]  # change "Date" if needed

png("BSSSM_PseudoResidualDiagnostics.png", width = 2400, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.8, 0))

col_main <- "#901a1E"
col_ref  <- "steelblue"

# (1,1) Time series
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

# (1,2) Histogram (trimmed 1–99%) + N(0,1)
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

# (2,1) QQ-plot
qqnorm(res_SSM_finite,
       main = "QQ-plot",
       col  = col_main,
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles",
       cex.main = 1.4)
qqline(res_SSM_finite, col = col_ref, lwd = 2)

# (2,2) ACF
acf_obj <- acf(res_SSM_finite, plot = FALSE)
plot(acf_obj,
     main = "",
     xlab = "Lag",
     ylab = "Autocorrelation",
     col  = col_main,
     lwd  = 2)
title("ACF", cex.main = 1.4)

mtext("BS-SSM Pseudo-Residuals",
      side  = 3,
      outer = TRUE,
      cex   = 1,
      font  = 2)

dev.off()
