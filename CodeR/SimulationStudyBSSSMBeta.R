# Continuous State-Space Simulation: BS–SSM_beta 

# ---- Path ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(ggplot2)
library(matrixStats)

# ---- Time step ----
dt <- 1/252

# ---- True parameters ----
n               <- 25000
mu0_true        <- 0.05
beta_mu_true    <- 0.20
sigma_true      <- 0.15
beta_sigma_true <- 0.60
phi_s_true      <- 0.98
sigma_s_true    <- 0.10
S0              <- 100

# ---- Simulation ----
simulate_SSM_BS_beta <- function(n, mu0_true, beta_mu_true,
                                 sigma_true, beta_sigma_true,
                                 phi_s_true, sigma_s_true,
                                 S0, dt = 1/252) {
  C      <- numeric(n)
  logS   <- numeric(n)
  r      <- numeric(n - 1)
  
  C[1]    <- 0.0
  logS[1] <- log(S0)
  
  for (i in 2:n) {
    C[i] <- phi_s_true * C[i - 1] + sigma_s_true * rnorm(1)
    
    mu_i    <- mu0_true + beta_mu_true * C[i]
    sigma_i <- sigma_true * exp(beta_sigma_true * C[i])
    
    incr     <- (mu_i - 0.5 * sigma_i^2) * dt + sigma_i * sqrt(dt) * rnorm(1)
    logS[i]  <- logS[i - 1] + incr
    r[i - 1] <- incr
  }
  
  S <- exp(logS)
  list(r = r, C = C, logS = logS, S = S)
}

# ---- Negative log-likelihood (Kitagawa grid; scaled forward recursion) ----
mllk_SSM_BS_beta <- function(theta.star, x, m, bm, dt = 1/252) {
  phi_state   <- plogis(theta.star[1])
  sigma_state <- exp(theta.star[2])
  mu0         <- theta.star[3]
  sigma       <- exp(theta.star[4])
  beta_mu     <- theta.star[5]
  beta_sigma  <- theta.star[6]
  
  b      <- seq(-bm, bm, length = m + 1)
  h      <- b[2] - b[1]
  bstar  <- 0.5 * (b[-1] + b[-(m + 1)])
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  sd0   <- sigma_state / sqrt(1 - phi_state^2)
  delta <- h * dnorm(bstar, mean = 0, sd = sd0)
  
  mu_vec    <- mu0 + beta_mu * bstar
  sigma_vec <- sigma * exp(beta_sigma * bstar)
  mean_vec  <- (mu_vec - 0.5 * sigma_vec^2) * dt
  sd_vec    <- sigma_vec * sqrt(dt)
  
  Tn <- length(x)
  
  dens1 <- dnorm(x[1], mean = mean_vec, sd = sd_vec)
  foo   <- delta * dens1
  l     <- log(sum(foo))
  phi   <- foo / sum(foo)
  
  for (t in 2:Tn) {
    denst <- dnorm(x[t], mean = mean_vec, sd = sd_vec)
    foo   <- as.numeric(phi %*% Gamma) * denst
    st    <- sum(foo)
    l     <- l + log(st)
    phi   <- foo / st
  }
  
  return(-l)
}

# ---- Initial guess (transformed scale) ----
theta.start_beta <- c(
  qlogis(phi_s_true),
  log(sigma_s_true),
  mu0_true,
  log(sigma_true),
  beta_mu_true,
  beta_sigma_true
)

# ---- Fit settings ----
num_simulations <- 1
m_grid          <- 200
bm_grid         <- 3

fitted_params_SSM_BS_beta <- matrix(NA_real_, nrow = num_simulations, ncol = length(theta.start_beta))

# ---- Simulate and fit ----
for (i in seq_len(num_simulations)) {
  set.seed(i)
  cat(sprintf("[%d/%d]\n", i, num_simulations))
  
  sim <- simulate_SSM_BS_beta(
    n               = n,
    mu0_true        = mu0_true,
    beta_mu_true    = beta_mu_true,
    sigma_true      = sigma_true,
    beta_sigma_true = beta_sigma_true,
    phi_s_true      = phi_s_true,
    sigma_s_true    = sigma_s_true,
    S0              = S0,
    dt              = dt
  )
  r <- sim$r
  
  est_i <- tryCatch({
    fit <- nlm(
      f = mllk_SSM_BS_beta,
      p = theta.start_beta,
      x = r,
      m = m_grid,
      bm = bm_grid,
      dt = dt,
      print.level = 0
    )
    if (!is.finite(fit$minimum) || any(!is.finite(fit$estimate))) {
      stop("non-finite estimate from nlm")
    }
    fit$estimate
  }, error = function(e) {
    message(sprintf("Seed %d failed: %s — skipping and recording NA.", i, conditionMessage(e)))
    rep(NA_real_, length(theta.start_beta))
  })
  
  fitted_params_SSM_BS_beta[i, ] <- est_i
}

# ---- Back-transform ----
fitted_params_original_SSM_BS_beta <- cbind(
  phi        = plogis(fitted_params_SSM_BS_beta[, 1]),
  sigZ       = exp(fitted_params_SSM_BS_beta[, 2]),
  mu0        = fitted_params_SSM_BS_beta[, 3],
  sigma      = exp(fitted_params_SSM_BS_beta[, 4]),
  beta_mu    = fitted_params_SSM_BS_beta[, 5],
  beta_sigma = fitted_params_SSM_BS_beta[, 6]
)

SSM_BS_beta_means <- colMeans(fitted_params_original_SSM_BS_beta, na.rm = TRUE)
print(SSM_BS_beta_means)

save(fitted_params_SSM_BS_beta, file = "fitted_params_BSSSM_beta.RData")

# ---- Summary: means + quantiles ----
q_probs  <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
q_labels <- paste0("q", sub("\\.?0+$","", format(q_probs*100, trim = TRUE)))

df_beta <- fitted_params_original_SSM_BS_beta
X_beta  <- if (is.matrix(df_beta)) {
  storage.mode(df_beta) <- "double"; df_beta
} else {
  as.matrix(df_beta[sapply(df_beta, is.numeric)])
}

if (is.null(colnames(X_beta))) colnames(X_beta) <- paste0("V", seq_len(ncol(X_beta)))

means_beta  <- colMeans(X_beta, na.rm = TRUE)
quants_beta <- matrixStats::colQuantiles(X_beta, probs = q_probs, na.rm = TRUE)
colnames(quants_beta) <- q_labels

SSM_BS_beta_summary <- data.frame(
  parameter   = colnames(X_beta),
  mean        = as.numeric(means_beta),
  quants_beta,
  check.names = FALSE,
  row.names   = NULL
)

SSM_BS_beta_summary[1:2]

# ---- Plot one simulated path ----
S_beta <- sim$S
t_years <- (0:(length(S_beta)-1)) * dt

op <- options(scipen = 999)

df_plot <- data.frame(Years = t_years, StockPrice = S_beta)

ggplot(df_plot, aes(x = Years, y = StockPrice)) +
  geom_line(linewidth = 0.9, color = "#901a1E") +
  labs(x = "Years", y = "Stock Price") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

options(op)
