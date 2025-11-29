###############################################################################
# BS-HMM: State-wise tables for all 12 models (μ, σ, μσ families, 2–5 states)
# - Uses your parameter indexing from the mllk definitions
# - SEs & 95% CIs for μ_cap and σ via Hessian + Jacobian (like make_param_table)
# - HAC SEs & CIs for q_i; μ_tot = μ_cap + q (no SE/CI)
# - Outputs LaTeX tables in your preferred format
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(sandwich)
  library(MASS)   # for ginv fallback if Hessian is near-singular
})

dt <- 1/252

#-------------------------------- Load data ----------------------------------#
load("sp500_clean_with_divs.RData")   # expects: spx_div with column 'DivLog'
load("sp500_clean.RData")             # expects: spx with column 'Return'
load("fitted_params_BSHMM_2019.RData")# expects: models_BSHMM_2019 list

# Ensure return vector matches dividend data length
spx   <- spx %>% slice_head(n = nrow(spx_div))
x_ret <- as.numeric(spx$Return)

###############################################################################
# Viterbi decoding (GBM returns)
###############################################################################
viterbi <- function(x, mu, sigma, Gamma, delta, N) {
  dt <- 1/252
  n  <- length(x)
  
  allprobs <- matrix(1, nrow = n, ncol = N)
  ind      <- which(!is.na(x))
  
  for (j in 1:N) {
    allprobs[ind, j] <- dnorm(
      x[ind],
      mean = (mu[j] - 0.5 * sigma[j]^2) * dt,
      sd   = sqrt(dt * sigma[j]^2)
    )
  }
  
  xi  <- matrix(0, nrow = n, ncol = N)
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  
  for (t in 2:n) {
    foo      <- apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ]  <- foo / sum(foo)
  }
  
  iv    <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ])
  }
  iv
}

###############################################################################
# Estimate mean + HAC SE for dividend yield in one state group
###############################################################################
estimate_q_hac <- function(q_vec) {
  mod <- lm(q_vec ~ 1)
  bw  <- sandwich::bwNeweyWest(mod)
  vc  <- sandwich::NeweyWest(mod, lag = bw, prewhite = FALSE, adjust = TRUE)
  
  q_hat_daily <- coef(mod)[1]
  se_daily    <- sqrt(diag(vc))[1]
  
  q_hat_annual <- 252 * q_hat_daily
  se_annual    <- 252 * se_daily
  
  list(q_hat_annual = q_hat_annual, se_annual = se_annual)
}

###############################################################################
# Extract μ_cap, σ, Γ, δ + SEs & 95% CIs consistent with your mllk + make_param_table
###############################################################################
get_bshmm_params <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
  family <- match.arg(family)
  theta  <- fit$estimate
  H      <- fit$hessian
  P      <- length(theta)
  
  if (is.null(H)) stop("Hessian not stored in fit object.")
  
  # Covariance of working parameters (μ, logσ, η) — with ginv fallback
  Cov_phi <- tryCatch(
    solve(H),
    error = function(e) MASS::ginv(H)
  )
  
  # Identify indices according to your mllk parameterisations
  if (family == "mu") {
    # θ = ( μ_1,...,μ_N, log σ, η )
    mu_idx        <- 1:N
    log_sigma_idx <- N + 1
    eta_idx       <- (N + 2):P
    
    mu_hat        <- theta[mu_idx]
    log_sigma     <- theta[log_sigma_idx]
    sigma_hat     <- rep(exp(log_sigma), N)
    
  } else if (family == "sigma") {
    # θ = ( μ, log σ_1,...,log σ_N, η )
    mu_idx        <- 1
    log_sigma_idx <- 2:(N + 1)
    eta_idx       <- (N + 2):P
    
    mu_hat    <- rep(theta[mu_idx], N)       # common μ repeated per state
    log_sigma <- theta[log_sigma_idx]
    sigma_hat <- exp(log_sigma)
    
  } else { # "mu_sigma"
    # θ = ( μ_1,...,μ_N, log σ_1,...,log σ_N, η )
    mu_idx        <- 1:N
    log_sigma_idx <- (N + 1):(2 * N)
    eta_idx       <- (2 * N + 1):P
    
    mu_hat    <- theta[mu_idx]
    log_sigma <- theta[log_sigma_idx]
    sigma_hat <- exp(log_sigma)
  }
  
  # Transition matrix Γ exactly as in your mllk code
  eta   <- theta[eta_idx]
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution δ (eigen of Γ^T, like in mllk)
  eig   <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  # Delta-method Jacobian: natural params (μ_1..μ_N, σ_1..σ_N)
  # rows = 2N natural params, cols = P working params
  J <- matrix(0, nrow = 2 * N, ncol = P)
  
  if (family == "mu") {
    # μ_i depends only on θ[mu_idx[i]]
    for (i in 1:N) {
      J[i, mu_idx[i]] <- 1
    }
    # σ_i (all identical) depend on log σ (single index)
    for (i in 1:N) {
      J[N + i, log_sigma_idx] <- sigma_hat[i]   # dσ_i / d logσ = σ
    }
    
  } else if (family == "sigma") {
    # μ_i (all identical) depend on θ[1] (common μ)
    for (i in 1:N) {
      J[i, mu_idx] <- 1
    }
    # σ_i depend on their own log σ_i
    for (i in 1:N) {
      J[N + i, log_sigma_idx[i]] <- sigma_hat[i]
    }
    
  } else { # mu_sigma
    # μ_i depends only on its own μ_i
    for (i in 1:N) {
      J[i, mu_idx[i]] <- 1
    }
    # σ_i depends only on its own log σ_i
    for (i in 1:N) {
      J[N + i, log_sigma_idx[i]] <- sigma_hat[i]
    }
  }
  
  # Covariance on natural scale (μ_i, σ_i)
  Cov_nat <- J %*% Cov_phi %*% t(J)
  diagCov <- diag(Cov_nat)
  diagCov[diagCov < 0] <- NA_real_  # guard numerical negatives
  
  se_vec    <- sqrt(diagCov)
  se_mu     <- se_vec[1:N]
  se_sigma  <- se_vec[(N + 1):(2 * N)]
  
  z <- 1.96
  mu_lower    <- mu_hat    - z * se_mu
  mu_upper    <- mu_hat    + z * se_mu
  sigma_lower <- sigma_hat - z * se_sigma
  sigma_upper <- sigma_hat + z * se_sigma
  
  list(
    N           = N,
    family      = family,
    mu_cap      = mu_hat,
    se_mu_cap   = se_mu,
    mu_lower    = mu_lower,
    mu_upper    = mu_upper,
    sigma       = sigma_hat,
    se_sigma    = se_sigma,
    sigma_lower = sigma_lower,
    sigma_upper = sigma_upper,
    Gamma       = Gamma,
    delta       = delta
  )
}

###############################################################################
# For one model: estimate state-wise q_i, then combine with μ_cap and σ
###############################################################################
get_statewise_params_with_divs <- function(fit, N, family,
                                           model_label,
                                           spx, spx_div, x_ret) {
  # 1) Extract μ_cap, σ, Γ, δ with SEs/CIs
  pars <- get_bshmm_params(fit, N, family)
  
  # 2) Viterbi path
  v_path <- viterbi(
    x     = x_ret,
    mu    = pars$mu_cap,
    sigma = pars$sigma,
    Gamma = pars$Gamma,
    delta = pars$delta,
    N     = N
  )
  
  # 3) Attach states to dividends and estimate q_i with HAC Newey–West
  spx_div_states <- spx_div %>%
    mutate(State = v_path) %>%
    filter(!is.na(DivLog), !is.na(State))
  
  q_df <- spx_div_states %>%
    group_by(State) %>%
    summarise(
      N_days = n(),
      q_daily = list(DivLog),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(est = list(estimate_q_hac(q_daily))) %>%
    mutate(
      q_hat_annual = est$q_hat_annual,
      se_q_annual  = est$se_annual
    ) %>%
    ungroup() %>%
    mutate(
      q_lower = q_hat_annual - 1.96 * se_q_annual,
      q_upper = q_hat_annual + 1.96 * se_q_annual
    ) %>%
    dplyr::select(State, N_days, q_hat_annual, se_q_annual, q_lower, q_upper) %>%
    dplyr::arrange(State)
  
  states <- 1:N
  
  mu_df <- tibble(
    State    = states,
    mu_cap   = pars$mu_cap,
    se_mu    = pars$se_mu_cap,
    mu_lower = pars$mu_lower,
    mu_upper = pars$mu_upper
  )
  
  sigma_df <- tibble(
    State        = states,
    sigma        = pars$sigma,
    se_sigma     = pars$se_sigma,
    sigma_lower  = pars$sigma_lower,
    sigma_upper  = pars$sigma_upper
  )
  
  comb <- mu_df %>%
    left_join(q_df,     by = "State") %>%
    left_join(sigma_df, by = "State") %>%
    mutate(
      Model  = model_label,
      mu_tot = mu_cap + q_hat_annual
    ) %>%
    dplyr::arrange(State)
  
  comb
}

###############################################################################
# Formatting helpers for LaTeX – robust to missing SE/CI (NA ⇒ dashes)
###############################################################################
format_ci <- function(est, lower, upper, digits = 4) {
  fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)
  
  if (is.na(lower) || is.na(upper)) {
    return(paste0(fmt(est), "\\; (\\textemdash, \\textemdash)"))
  }
  paste0(
    fmt(est), "\\; (",
    fmt(lower), ", ",
    fmt(upper), ")"
  )
}

format_se <- function(se, digits = 5) {
  if (!is.finite(se)) {
    return("\\textemdash")
  }
  sprintf(paste0("%.", digits, "f"), se)
}

# Build rows: μ_cap for all states, then q, then μ_tot, then σ
build_param_rows <- function(state_df) {
  rows <- list()
  
  # 1) μ_cap,i rows
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{\\mu}_{\\text{cap},%d}$", s),
      Est_CI    = format_ci(r$mu_cap, r$mu_lower, r$mu_upper),
      SE        = format_se(r$se_mu)
    )
  }
  
  # 2) q_i rows
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{q}_{%d}$", s),
      Est_CI    = format_ci(r$q_hat_annual, r$q_lower, r$q_upper),
      SE        = format_se(r$se_q_annual)
    )
  }
  
  # 3) μ_tot,i rows (no SE / CI by design)
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{\\mu}_{\\text{tot},%d}$", s),
      Est_CI    = paste0(
        sprintf("%.4f", r$mu_tot),
        "\\; (\\textemdash, \\textemdash)"
      ),
      SE        = "\\textemdash"
    )
  }
  
  # 4) σ_i rows
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{\\sigma}_{%d}$", s),
      Est_CI    = format_ci(r$sigma, r$sigma_lower, r$sigma_upper),
      SE        = format_se(r$se_sigma)
    )
  }
  
  dplyr::bind_rows(rows)
}

build_latex_table <- function(param_rows, caption, label) {
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\begin{tabular}{lcc}",
    "\\toprule",
    "Parameter & Estimate (95\\% CI) & Std.\\ Error \\\\",
    "\\midrule"
  )
  
  for (i in seq_len(nrow(param_rows))) {
    lines <- c(
      lines,
      sprintf("%s & %s & %s \\\\",
              param_rows$Parameter[i],
              param_rows$Est_CI[i],
              param_rows$SE[i])
    )
  }
  
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\end{table}"
  )
  
  paste(lines, collapse = "\n")
}

###############################################################################
# Model list: 12 models in models_BSHMM_2019
###############################################################################
model_specs <- tibble::tribble(
  ~name,           ~N, ~family,
  "mod2_mu",        2, "mu",
  "mod3_mu",        3, "mu",
  "mod4_mu",        4, "mu",
  "mod5_mu",        5, "mu",
  "mod2_sigma",     2, "sigma",
  "mod3_sigma",     3, "sigma",
  "mod4_sigma",     4, "sigma",
  "mod5_sigma",     5, "sigma",
  "mod2_mu_sigma",  2, "mu_sigma",
  "mod3_mu_sigma",  3, "mu_sigma",
  "mod4_mu_sigma",  4, "mu_sigma",
  "mod5_mu_sigma",  5, "mu_sigma"
)

family_desc <- list(
  mu       = "state-dependent drift $\\mu_{\\text{cap},i}$ with common volatility",
  sigma    = "state-dependent volatility $\\sigma_i$ with common drift",
  mu_sigma = "state-dependent drift $\\mu_{\\text{cap},i}$ and volatility $\\sigma_i$"
)

latex_tables <- list()

###############################################################################
# Main loop: build tables for all 12 models
###############################################################################
for (k in seq_len(nrow(model_specs))) {
  spec <- model_specs[k, ]
  fit  <- models_BSHMM_2019[[ spec$name ]]
  
  message("Processing model: ", spec$name)
  
  # 1) State-wise μ_cap, σ, q, μ_tot
  state_df <- get_statewise_params_with_divs(
    fit         = fit,
    N           = spec$N,
    family      = spec$family,
    model_label = spec$name,
    spx         = spx,
    spx_div     = spx_div,
    x_ret       = x_ret
  )
  
  # 2) Build rows in your preferred order
  param_rows <- build_param_rows(state_df)
  
  # 3) Caption + label
  cap <- sprintf(
    "%d-state BS-HMM with %s. State-dependent parameter estimates (annualised) with 95\\%% confidence intervals based on the inverse Hessian of the minimised log-likelihood on the working scale, transformed to the natural scale by the delta method. The capital-gains drift $\\mu_{\\text{cap},i}$ and dividend yield $q_i$ combine to give the total-return drift $\\mu_{\\text{tot},i} = \\mu_{\\text{cap},i} + q_i$ in each state. Standard errors and confidence intervals for $\\mu_{\\text{tot},i}$ are not reported and are indicated by dashes.",
    spec$N,
    family_desc[[ spec$family ]]
  )
  
  lab <- sprintf("tab:BSHMM_%s", spec$name)
  
  latex_tables[[ spec$name ]] <- build_latex_table(param_rows, caption = cap, label = lab)
}

###############################################################################
# Print all tables to console (in LaTeX format)
###############################################################################
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}

###############################################################################
# Optional: write all tables to a single .tex file
###############################################################################
sink("BSHMM_state_tables_all_12_models.tex")
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}
sink()
