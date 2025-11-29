## ---- AR(1) MLE via closed-form (conditional likelihood) ----------------
## Model: x_t = ρ x_{t-1} + ε_t,  ε_t ~ N(0, σ^2)
## MLEs: ρ̂ = Σ x_t x_{t-1} / Σ x_{t-1}^2 ,   σ̂^2 = (1/N) Σ (x_t - ρ̂ x_{t-1})^2
## where sums are t = 2..T and N = T-1.

ar1_mle <- function(x, use_full_T = FALSE) {
  x <- as.numeric(stats::na.omit(x))
  if (length(x) < 2L) stop("Need at least 2 observations.")
  
  xt   <- x[-1L]
  xtm1 <- x[-length(x)]
  N    <- length(xt)            # = T - 1 pairs
  
  rho_hat    <- sum(xt * xtm1) / sum(xtm1^2)
  denom_T    <- if (use_full_T) length(x) else N
  sigma2_hat <- sum((xt - rho_hat * xtm1)^2) / denom_T
  
  # Likelihood (ignoring constants), matching your display: L ∝ (σ̂^2)^(-denom_T/2)
  L_hat      <- sigma2_hat^(-denom_T/2)
  logLik_hat <- -(denom_T/2) * (1 + log(2*pi*sigma2_hat))
  
  list(T_total = length(x),
       N_pairs = N,
       rho_hat = rho_hat,
       sigma2_hat = sigma2_hat,
       L_hat = L_hat,
       logLik_ign_const = logLik_hat)
}

## ---- Run on your two series --------------------------------------------
res_prices  <- ar1_mle(spx$AdjClose)     # prices (often near unit root)
res_returns <- ar1_mle(spx$Return)       # log returns (typically stationary-ish)

res_prices
##############
res_returns
