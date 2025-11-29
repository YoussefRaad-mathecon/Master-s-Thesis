###############################################################################
# BS-HMM: State-wise tables for all 12 models (μ, σ, μσ families, 2–5 states)
# - Uses your parameter indexing from the mllk definitions
# - Includes μ_cap, q, μ_tot = μ_cap + q, and σ (annualised where relevant)
# - Handles singular / ill-conditioned Hessians by skipping SE/CI (dashes)
# - Prints LaTeX tables for all models
###############################################################################

library(tidyverse)
library(sandwich)

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
# Extract parameters from nlm fit in a way consistent with YOUR mllk code,
# with robust handling of singular / bad Hessians
###############################################################################
get_bshmm_params <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
  family <- match.arg(family)
  th     <- fit$estimate
  H      <- fit$hessian
  
  if (is.null(H)) {
    warning("Hessian not stored in fit object; SEs and CIs will be NA.")
    se <- rep(NA_real_, length(th))
  } else {
    V <- tryCatch(
      solve(H),
      error = function(e) {
        warning("Hessian is singular or not invertible; SEs and CIs will be NA.")
        NULL
      }
    )
    if (is.null(V)) {
      se <- rep(NA_real_, length(th))
    } else {
      diagV <- diag(V)
      # Negative variances → treat as NA (no SE)
      diagV[diagV < 0] <- NA_real_
      se <- sqrt(diagV)
    }
  }
  
  if (family == "mu") {
    # θ = ( μ_1,...,μ_N, log σ, η )
    mu_idx        <- 1:N
    log_sigma_idx <- N + 1
    eta_idx       <- (N + 2):length(th)
    
    mu_raw        <- th[mu_idx]
    se_mu_raw     <- se[mu_idx]
    
    log_sigma_raw <- th[log_sigma_idx]
    se_log_sigma  <- se[log_sigma_idx]
    
    log_sigma     <- rep(log_sigma_raw, N)
    se_log_sigma  <- rep(se_log_sigma, N)
    
  } else if (family == "sigma") {
    # θ = ( μ, log σ_1,...,log σ_N, η )
    mu_idx        <- 1
    log_sigma_idx <- 2:(N + 1)
    eta_idx       <- (N + 2):length(th)
    
    mu_raw        <- rep(th[mu_idx], N)
    se_mu_raw     <- rep(se[mu_idx], N)
    
    log_sigma     <- th[log_sigma_idx]
    se_log_sigma  <- se[log_sigma_idx]
    
  } else { # "mu_sigma"
    # θ = ( μ_1,...,μ_N, log σ_1,...,log σ_N, η )
    mu_idx        <- 1:N
    log_sigma_idx <- (N + 1):(2 * N)
    eta_idx       <- (2 * N + 1):length(th)
    
    mu_raw        <- th[mu_idx]
    se_mu_raw     <- se[mu_idx]
    
    log_sigma     <- th[log_sigma_idx]
    se_log_sigma  <- se[log_sigma_idx]
  }
  
  eta <- th[eta_idx]
  
  # σ and its SEs via delta method (might be NA)
  sigma        <- exp(log_sigma)
  se_sigma     <- sigma * se_log_sigma
  
  # 95% CIs (become NA if SEs are NA)
  mu_lower     <- mu_raw   - 1.96 * se_mu_raw
  mu_upper     <- mu_raw   + 1.96 * se_mu_raw
  
  sigma_lower  <- exp(log_sigma - 1.96 * se_log_sigma)
  sigma_upper  <- exp(log_sigma + 1.96 * se_log_sigma)
  
  # Transition matrix Γ exactly as in your code
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  
  # Stationary distribution δ via eigen (same as in your mllk)
  eig   <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  list(
    N           = N,
    family      = family,
    mu_cap      = mu_raw,
    se_mu_cap   = se_mu_raw,
    mu_lower    = mu_lower,
    mu_upper    = mu_upper,
    sigma       = sigma,
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
    select(State, N_days, q_hat_annual, se_q_annual, q_lower, q_upper) %>%
    arrange(State)
  
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
    arrange(State)
  
  comb
}

###############################################################################
# Formatting helpers for LaTeX – robust to missing SE/CI (NA ⇒ dashes)
###############################################################################
format_ci <- function(est, lower, upper, digits = 4) {
  fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)
  
  if (any(is.na(c(lower, upper)))) {
    # Point estimate available, but CI not → show (—, —)
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
  
  # 2) q_i rows (HAC SEs should be fine)
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
    "%d-state BS-HMM with %s. State-dependent parameter estimates (annualised) with 95\\%% confidence intervals based on the inverse Hessian of the minimised log-likelihood where available. The capital-gains drift $\\mu_{\\text{cap},i}$ and dividend yield $q_i$ combine to give the total-return drift $\\mu_{\\text{tot},i} = \\mu_{\\text{cap},i} + q_i$ in each state. Standard errors and confidence intervals for $\\mu_{\\text{tot},i}$ are not reported and are indicated by dashes. Where the Hessian is singular or ill-conditioned, standard errors and confidence intervals for $\\mu_{\\text{cap},i}$ and $\\sigma_i$ are omitted and indicated by dashes.",
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
