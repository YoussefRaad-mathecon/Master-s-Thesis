# BS–HMM: State-wise tables + Markov parameters for 12 models

# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(tidyverse)
library(sandwich)
library(MASS)

# ---- Settings ----
dt <- 1/252

# ---- Data ----
load("sp500_clean_with_divs.RData")    # spx_div with 'DivLog'
load("sp500_clean.RData")              # spx with 'Return'
load("fitted_params_BSHMM_2019.RData") # models_BSHMM_2019

spx   <- spx %>% slice_head(n = nrow(spx_div))
x_ret <- as.numeric(spx$Return)

# ---- Viterbi decoding (GBM returns) ----
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

# ---- HAC mean + SE for dividend yield ----
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

# ---- Core BS–HMM parameters + delta-method SEs on natural scale ----
get_bshmm_params <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
  family <- match.arg(family)
  theta  <- fit$estimate
  H      <- fit$hessian
  P      <- length(theta)
  
  if (is.null(H)) stop("Hessian not stored in fit object.")
  
  Cov_phi <- tryCatch(
    solve(H),
    error = function(e) MASS::ginv(H)
  )
  
  if (family == "mu") {
    mu_idx        <- 1:N
    log_sigma_idx <- N + 1
    eta_idx       <- (N + 2):P
    
    mu_hat        <- theta[mu_idx]
    log_sigma     <- theta[log_sigma_idx]
    sigma_hat     <- rep(exp(log_sigma), N)
    
  } else if (family == "sigma") {
    mu_idx        <- 1
    log_sigma_idx <- 2:(N + 1)
    eta_idx       <- (N + 2):P
    
    mu_hat    <- rep(theta[mu_idx], N)
    log_sigma <- theta[log_sigma_idx]
    sigma_hat <- exp(log_sigma)
    
  } else {
    mu_idx        <- 1:N
    log_sigma_idx <- (N + 1):(2 * N)
    eta_idx       <- (2 * N + 1):P
    
    mu_hat    <- theta[mu_idx]
    log_sigma <- theta[log_sigma_idx]
    sigma_hat <- exp(log_sigma)
  }
  
  eta   <- theta[eta_idx]
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  
  eig   <- eigen(t(Gamma))
  delta <- Re(eig$vectors[, 1])
  delta <- delta / sum(delta)
  
  J <- matrix(0, nrow = 2 * N, ncol = P)
  
  if (family == "mu") {
    for (i in 1:N) {
      J[i, mu_idx[i]] <- 1
    }
    for (i in 1:N) {
      J[N + i, log_sigma_idx] <- sigma_hat[i]
    }
    
  } else if (family == "sigma") {
    for (i in 1:N) {
      J[i, mu_idx] <- 1
    }
    for (i in 1:N) {
      J[N + i, log_sigma_idx[i]] <- sigma_hat[i]
    }
    
  } else {
    for (i in 1:N) {
      J[i, mu_idx[i]] <- 1
    }
    for (i in 1:N) {
      J[N + i, log_sigma_idx[i]] <- sigma_hat[i]
    }
  }
  
  Cov_nat <- J %*% Cov_phi %*% t(J)
  diagCov <- diag(Cov_nat)
  diagCov[diagCov < 0] <- NA_real_
  
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
    delta       = delta,
    theta       = theta,
    Cov_phi     = Cov_phi,
    eta_idx     = eta_idx
  )
}

# ---- Markov parameters + delta-method SEs (numeric Jacobian) ----
get_markov_params <- function(pars) {
  N       <- pars$N
  theta   <- pars$theta
  Cov_phi <- pars$Cov_phi
  eta_idx <- pars$eta_idx
  
  eta_hat <- theta[eta_idx]
  Cov_eta <- Cov_phi[eta_idx, eta_idx, drop = FALSE]
  
  markov_natural <- function(eta_vec) {
    Gamma <- diag(N)
    Gamma[!Gamma] <- exp(eta_vec)
    Gamma <- Gamma / rowSums(Gamma)
    
    eig   <- eigen(t(Gamma))
    delta <- Re(eig$vectors[, 1])
    delta <- delta / sum(delta)
    
    c(as.vector(Gamma), delta)
  }
  
  psi0 <- markov_natural(eta_hat)
  Q    <- length(psi0)
  Peta <- length(eta_hat)
  
  J_eta <- matrix(0, nrow = Q, ncol = Peta)
  for (k in 1:Peta) {
    step        <- 1e-6 * max(1, abs(eta_hat[k]))
    eta_step    <- eta_hat
    eta_step[k] <- eta_step[k] + step
    psi_step    <- markov_natural(eta_step)
    J_eta[, k]  <- (psi_step - psi0) / step
  }
  
  Cov_markov <- J_eta %*% Cov_eta %*% t(J_eta)
  diagCov    <- diag(Cov_markov)
  diagCov[diagCov < 0] <- NA_real_
  se_vec     <- sqrt(diagCov)
  
  n_gamma      <- N^2
  se_gamma_vec <- se_vec[1:n_gamma]
  se_delta     <- se_vec[(n_gamma + 1):(n_gamma + N)]
  
  Gamma_hat <- matrix(psi0[1:n_gamma], nrow = N, ncol = N)
  se_Gamma  <- matrix(se_gamma_vec,    nrow = N, ncol = N)
  delta_hat <- psi0[(n_gamma + 1):(n_gamma + N)]
  
  list(
    Gamma    = Gamma_hat,
    se_Gamma = se_Gamma,
    delta    = delta_hat,
    se_delta = se_delta
  )
}

# ---- State-wise q_i + join with mu_cap and sigma ----
get_statewise_params_with_divs <- function(fit, N, family,
                                           model_label,
                                           spx, spx_div, x_ret,
                                           pars = NULL) {
  if (is.null(pars)) {
    pars <- get_bshmm_params(fit, N, family)
  }
  
  v_path <- viterbi(
    x     = x_ret,
    mu    = pars$mu_cap,
    sigma = pars$sigma,
    Gamma = pars$Gamma,
    delta = pars$delta,
    N     = N
  )
  
  spx_div_states <- spx_div %>%
    mutate(State = v_path) %>%
    filter(!is.na(DivLog), !is.na(State))
  
  q_df <- spx_div_states %>%
    group_by(State) %>%
    summarise(
      N_days  = n(),
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

# ---- LaTeX formatting helpers ----
format_ci <- function(est, lower, upper, digits = 4) {
  fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)
  
  if (is.na(lower) || is.na(upper)) {
    return(paste0(fmt(est), "\\; (\\textemdash, \\textemdash)"))
  }
  paste0(fmt(est), "\\; (", fmt(lower), ", ", fmt(upper), ")")
}

format_se <- function(se, digits = 5) {
  if (!is.finite(se)) {
    return("\\textemdash")
  }
  sprintf(paste0("%.", digits, "f"), se)
}

# ---- Build LaTeX table rows ----
build_param_rows <- function(state_df) {
  rows <- list()
  
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{\\mu}_{\\text{cap},%d}$", s),
      Est_CI    = format_ci(r$mu_cap, r$mu_lower, r$mu_upper),
      SE        = format_se(r$se_mu)
    )
  }
  
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{q}_{%d}$", s),
      Est_CI    = format_ci(r$q_hat_annual, r$q_lower, r$q_upper),
      SE        = format_se(r$se_q_annual)
    )
  }
  
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble::tibble(
      Parameter = sprintf("$\\widehat{\\mu}_{\\text{tot},%d}$", s),
      Est_CI    = paste0(sprintf("%.4f", r$mu_tot), "\\; (\\textemdash, \\textemdash)"),
      SE        = "\\textemdash"
    )
  }
  
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

# ---- LaTeX blocks for Gamma and delta ----
matrix_to_pmatrix <- function(M, digits = 4) {
  rows <- apply(M, 1, function(row) {
    paste(sprintf(paste0("%.", digits, "f"), row), collapse = " & ")
  })
  paste0("\\begin{pmatrix}\n",
         paste(rows, collapse = " \\\\ \n"),
         "\n\\end{pmatrix}")
}

vector_to_parens <- function(v, digits = 4) {
  paste0("\\left(",
         paste(sprintf(paste0("%.", digits, "f"), v), collapse = ", "),
         "\\right)")
}

build_markov_tex_block <- function(markov) {
  Gamma_hat_tex <- matrix_to_pmatrix(markov$Gamma,    digits = 4)
  se_Gamma_tex  <- matrix_to_pmatrix(markov$se_Gamma, digits = 4)
  delta_tex     <- vector_to_parens(markov$delta,     digits = 4)
  se_delta_tex  <- vector_to_parens(markov$se_delta,  digits = 4)
  
  c(
    sprintf("\\multicolumn{3}{l}{$\\widehat{\\Gamma} = %s,\\quad \\widehat{\\mathrm{se}}(\\widehat{\\Gamma}) = %s$}\\\\",
            Gamma_hat_tex, se_Gamma_tex),
    sprintf("\\multicolumn{3}{l}{$\\widehat{\\delta} = %s,\\quad \\widehat{\\mathrm{se}}(\\widehat{\\delta}) = %s$}\\\\",
            delta_tex, se_delta_tex)
  )
}

# ---- LaTeX table builder ----
build_latex_table <- function(param_rows, caption, label, markov_tex = NULL) {
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
  
  if (!is.null(markov_tex)) {
    lines <- c(lines, "\\midrule", markov_tex)
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

# ---- Models ----
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

# ---- Main loop ----
for (k in seq_len(nrow(model_specs))) {
  spec <- model_specs[k, ]
  fit  <- models_BSHMM_2019[[ spec$name ]]
  
  message("Processing model: ", spec$name)
  
  pars <- get_bshmm_params(fit, N = spec$N, family = spec$family)
  
  state_df <- get_statewise_params_with_divs(
    fit         = fit,
    N           = spec$N,
    family      = spec$family,
    model_label = spec$name,
    spx         = spx,
    spx_div     = spx_div,
    x_ret       = x_ret,
    pars        = pars
  )
  
  param_rows <- build_param_rows(state_df)
  
  markov     <- get_markov_params(pars)
  markov_tex <- build_markov_tex_block(markov)
  
  cap <- sprintf(
    "%d-state BS-HMM with %s. State-dependent parameter estimates (annualised) with 95\\%% confidence intervals based on the inverse Hessian of the minimised log-likelihood on the working scale, transformed to the natural scale by the delta method. The capital-gains drift $\\mu_{\\text{cap},i}$ and dividend yield $q_i$ combine to give the total-return drift $\\mu_{\\text{tot},i} = \\mu_{\\text{cap},i} + q_i$ in each state. The bottom block reports the estimated transition matrix $\\boldsymbol{\\Gamma}$ and stationary distribution $\\boldsymbol{\\delta}$ together with their asymptotic standard errors.",
    spec$N,
    family_desc[[ spec$family ]]
  )
  
  lab <- sprintf("tab:BSHMM_%s", spec$name)
  
  latex_tables[[ spec$name ]] <- build_latex_table(
    param_rows,
    caption    = cap,
    label      = lab,
    markov_tex = markov_tex
  )
}

# ---- Print all tables ----
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}

# ---- Write all tables to file ----
sink("BSHMM_state_and_markov_tables_all_12_models.tex")
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}
sink()
