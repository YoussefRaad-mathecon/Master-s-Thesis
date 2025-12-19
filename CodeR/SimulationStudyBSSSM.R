# Continuous State-Space Simulation: BS–SSM (GBM)

# ---- Path ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(ggplot2)
library(matrixStats)

# ---- Time step ----
dt <- 1/252

# ---- True parameters ----
n            <- 25000
mu_true      <- 0.05
sigma_true   <- 0.15
phi_s_true   <- 0.98
sigma_s_true <- 0.10
S0           <- 100

# ---- Simulation: latent AR(1) + state-dependent drift/vol ----
simulate_SSM_BS <- function(n, mu_true, sigma_true, phi_s_true, sigma_s_true, S0, dt = 1/252) {
  Z      <- numeric(n)
  logS   <- numeric(n)
  r      <- numeric(n - 1)
  
  Z[1]    <- 0.0
  logS[1] <- log(S0)
  
  for (i in 2:n) {
    Z[i] <- phi_s_true * Z[i - 1] + sigma_s_true * rnorm(1)
    
    mu_i    <- mu_true    * exp(Z[i])
    sigma_i <- sigma_true * exp(Z[i])
    
    incr     <- (mu_i - 0.5 * sigma_i^2) * dt + sigma_i * sqrt(dt) * rnorm(1)
    logS[i]  <- logS[i - 1] + incr
    r[i - 1] <- incr
  }
  
  S <- exp(logS)
  list(r = r, Z = Z, logS = logS, S = S)
}

# ---- Negative log-likelihood (Kitagawa grid; scaled forward recursion) ----
mllk_SSM_BS <- function(theta.star, x, m, bm, dt = 1/252) {
  phi_state   <- plogis(theta.star[1])
  sigma_state <- exp(theta.star[2])
  mu          <- theta.star[3]
  sigma       <- exp(theta.star[4])
  
  b      <- seq(-bm, bm, length = m + 1)
  h      <- b[2] - b[1]
  bstar  <- 0.5 * (b[-1] + b[-(m + 1)])
  
  Gamma <- matrix(0, m, m)
  for (i in 1:m) {
    Gamma[i, ] <- h * dnorm(bstar, mean = phi_state * bstar[i], sd = sigma_state)
  }
  
  delta <- h * dnorm(bstar, mean = 0, sd = sigma_state / sqrt(1 - phi_state^2))
  
  mean_vec <- (mu * exp(bstar) - 0.5 * (sigma * exp(bstar))^2) * dt
  sd_vec   <- (sigma * exp(bstar)) * sqrt(dt)
  
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
theta.start <- c(
  qlogis(phi_s_true),
  log(sigma_s_true),
  mu_true,
  log(sigma_true)
)

# ---- Fit settings ----
num_simulations <- 1
m_grid          <- 200
bm_grid         <- 3

fitted_params_SSM_BS <- matrix(NA_real_, nrow = num_simulations, ncol = length(theta.start))

# ---- Simulate and fit ----
for (i in seq_len(num_simulations)) {
  set.seed(i)
  cat(sprintf("[%d/%d]\n", i, num_simulations))
  
  sim <- simulate_SSM_BS(n, mu_true, sigma_true, phi_s_true, sigma_s_true, S0, dt)
  r   <- sim$r
  
  est_i <- tryCatch({
    fit <- nlm(
      f = mllk_SSM_BS,
      p = theta.start,
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
    rep(NA_real_, length(theta.start))
  })
  
  fitted_params_SSM_BS[i, ] <- est_i
}

# ---- Back-transform ----
fitted_params_original_SSM_BS <- cbind(
  phi   = plogis(fitted_params_SSM_BS[, 1]),
  sigZ  = exp(fitted_params_SSM_BS[, 2]),
  mu    = fitted_params_SSM_BS[, 3],
  sigma = exp(fitted_params_SSM_BS[, 4])
)

SSM_BS_means <- colMeans(fitted_params_original_SSM_BS, na.rm = TRUE)
print(SSM_BS_means)

save(fitted_params_SSM_BS, file = "fitted_params_BSSSM.RData")

# ---- Summary: means + quantiles ----
q_probs  <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
q_labels <- paste0("q", sub("\\.?0+$","", format(q_probs*100, trim = TRUE)))

df <- fitted_params_original_SSM_BS
X  <- if (is.matrix(df)) {
  storage.mode(df) <- "double"; df
} else {
  as.matrix(df[sapply(df, is.numeric)])
}

if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))

means  <- colMeans(X, na.rm = TRUE)
quants <- matrixStats::colQuantiles(X, probs = q_probs, na.rm = TRUE)
colnames(quants) <- q_labels

SSM_BS_summary <- data.frame(
  parameter   = colnames(X),
  mean        = as.numeric(means),
  quants,
  check.names = FALSE,
  row.names   = NULL
)

SSM_BS_summary[1:2]

# ---- Plot helper ----
plot_price_path <- function(S, dt, title = NULL) {
  t_years <- (0:(length(S) - 1)) * dt
  
  df_plot <- data.frame(
    Years      = t_years,
    StockPrice = S
  )
  
  ggplot(df_plot, aes(x = Years, y = StockPrice)) +
    geom_line(linewidth = 0.9, color = "#901a1E") +
    labs(
      x = "Years",
      y = "Stock Price",
      title = title
    ) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 14, face = "plain", hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "plain"),
      axis.title.y = element_text(size = 14, face = "plain"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
}

# ---- Plot one simulated path ----
op <- options(scipen = 999)

S <- sim$S
p_bsssm <- plot_price_path(S, dt, title = NULL)
print(p_bsssm)

options(op)
