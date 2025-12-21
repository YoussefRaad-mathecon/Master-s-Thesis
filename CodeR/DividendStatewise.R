# BS-HMM: State-wise tables for all 12 models (μ, σ, μσ families, 2–5 states)


# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(tidyverse)
library(sandwich)
library(purrr)
library(tidyr)

# ---- Setup ----
set.seed(112233)
dt <- 1/252

# ---- Load data ----
load("sp500_clean_with_divs.RData")    # expects: spx_div with column 'DivLog'
load("sp500_clean.RData")              # expects: spx with column 'Return'
load("fitted_params_BSHMM_2019.RData") # expects: models_BSHMM_2019 list

stopifnot(exists("spx_div"), exists("spx"), exists("models_BSHMM_2019"))
stopifnot("DivLog" %in% names(spx_div), "Return" %in% names(spx))

# Ensure return vector matches dividend data length
spx   <- spx %>% dplyr::slice_head(n = nrow(spx_div))
x_ret <- as.numeric(spx$Return)

# ---- Viterbi decoding ----
viterbi <- function(x, mu, sigma, Gamma, delta, N, dt = 1/252) {
  n <- length(x)

  allprobs <- matrix(1, nrow = n, ncol = N)
  ind      <- which(!is.na(x))

  for (j in 1:N) {
    allprobs[ind, j] <- dnorm(
      x[ind],
      mean = (mu[j] - 0.5 * sigma[j]^2) * dt,
      sd   = sigma[j] * sqrt(dt)
    )
  }

  xi  <- matrix(0, nrow = n, ncol = N)
  foo <- delta * allprobs[1, ]
  s0  <- sum(foo)
  if (!is.finite(s0) || s0 <= 0) stop("Invalid initial scaling in Viterbi.")
  xi[1, ] <- foo / s0

  for (t in 2:n) {
    trans_weighted <- sweep(Gamma, 1, xi[t - 1, ], "*")
    foo            <- apply(trans_weighted, 2, max) * allprobs[t, ]
    st             <- sum(foo)
    if (!is.finite(st) || st <= 0) stop("Invalid scaling in Viterbi recursion at t = ", t)
    xi[t, ]        <- foo / st
  }

  iv    <- integer(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ])
  }
  iv
}

# ---- Dividend yield (robust HAC per state) ----
estimate_q_hac <- function(q_vec, min_n_hac = 30L, iid_fallback = FALSE) {
  q_vec <- as.numeric(q_vec)
  q_vec <- q_vec[is.finite(q_vec)]
  n     <- length(q_vec)

  if (n == 0L) {
    return(list(q_hat_annual = NA_real_, se_annual = NA_real_, lag = NA_integer_, n = 0L))
  }

  # Point estimate (annualised)
  q_hat_daily  <- mean(q_vec)
  q_hat_annual <- 252 * q_hat_daily

  # Too few obs for HAC: return NA SE (or optional IID fallback)
  if (n < min_n_hac) {
    se_daily <- if (iid_fallback && n >= 2L) stats::sd(q_vec) / sqrt(n) else NA_real_
    return(list(q_hat_annual = q_hat_annual, se_annual = 252 * se_daily, lag = 0L, n = n))
  }

  mod <- lm(q_vec ~ 1)

  bw <- tryCatch(sandwich::bwNeweyWest(mod), error = function(e) NA_real_)
  lag <- if (is.finite(bw)) as.integer(round(bw)) else 0L
  lag <- max(0L, min(lag, n - 1L))

  vc <- tryCatch(
    sandwich::NeweyWest(mod, lag = lag, prewhite = FALSE, adjust = TRUE),
    error = function(e) NULL
  )

  if (is.null(vc) || !is.finite(vc[1, 1]) || vc[1, 1] < 0) {
    se_daily <- if (iid_fallback && n >= 2L) stats::sd(q_vec) / sqrt(n) else NA_real_
    return(list(q_hat_annual = q_hat_annual, se_annual = 252 * se_daily, lag = 0L, n = n))
  }

  se_daily <- sqrt(vc[1, 1])
  list(q_hat_annual = q_hat_annual, se_annual = 252 * se_daily, lag = lag, n = n)
}

# ---- Extract parameters from nlm fit (consistent with your mllk code) ----
get_bshmm_params <- function(fit, N, family = c("mu", "sigma", "mu_sigma")) {
  family <- match.arg(family)
  th     <- fit$estimate
  H      <- fit$hessian

  if (is.null(th)) stop("Fit object has no $estimate.")

  # SEs on working scale (inverse Hessian), if possible
  if (is.null(H)) {
    warning("Hessian not stored; SEs and CIs set to NA.")
    se <- rep(NA_real_, length(th))
  } else {
    V <- tryCatch(
      solve(H),
      error = function(e) NULL
    )
    if (is.null(V)) {
      warning("Hessian singular/non-invertible; SEs and CIs set to NA.")
      se <- rep(NA_real_, length(th))
    } else {
      dv <- diag(V)
      dv[dv < 0] <- NA_real_
      se <- sqrt(dv)
    }
  }

  if (family == "mu") {
    mu_idx        <- 1:N
    log_sigma_idx <- N + 1
    eta_idx       <- (N + 2):length(th)

    mu_raw       <- th[mu_idx]
    se_mu_raw    <- se[mu_idx]
    log_sigma    <- rep(th[log_sigma_idx], N)
    se_log_sigma <- rep(se[log_sigma_idx], N)

  } else if (family == "sigma") {
    mu_idx        <- 1
    log_sigma_idx <- 2:(N + 1)
    eta_idx       <- (N + 2):length(th)

    mu_raw       <- rep(th[mu_idx], N)
    se_mu_raw    <- rep(se[mu_idx], N)
    log_sigma    <- th[log_sigma_idx]
    se_log_sigma <- se[log_sigma_idx]

  } else { # mu_sigma
    mu_idx        <- 1:N
    log_sigma_idx <- (N + 1):(2 * N)
    eta_idx       <- (2 * N + 1):length(th)

    mu_raw       <- th[mu_idx]
    se_mu_raw    <- se[mu_idx]
    log_sigma    <- th[log_sigma_idx]
    se_log_sigma <- se[log_sigma_idx]
  }

  eta <- th[eta_idx]

  sigma    <- exp(log_sigma)
  se_sigma <- sigma * se_log_sigma

  z <- 1.96
  mu_lower <- mu_raw - z * se_mu_raw
  mu_upper <- mu_raw + z * se_mu_raw

  sigma_lower <- exp(log_sigma - z * se_log_sigma)
  sigma_upper <- exp(log_sigma + z * se_log_sigma)

  # Transition matrix from eta (off-diagonals)
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)

  # Stationary distribution via eigenvector of t(Gamma)
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

# ---- State-wise q_i + combine with mu_cap and sigma ----
get_statewise_params_with_divs <- function(fit, N, family, model_label, spx_div, x_ret, dt = 1/252) {
  pars <- get_bshmm_params(fit, N, family)

  v_path <- viterbi(
    x     = x_ret,
    mu    = pars$mu_cap,
    sigma = pars$sigma,
    Gamma = pars$Gamma,
    delta = pars$delta,
    N     = N,
    dt    = dt
  )

  spx_div_states <- spx_div %>%
    mutate(State = v_path) %>%
    filter(!is.na(DivLog), !is.na(State))

  q_df <- spx_div_states %>%
    group_by(State) %>%
    summarise(
      N_days      = n(),
      q_daily_vec = list(as.numeric(DivLog)),
      .groups     = "drop"
    ) %>%
    mutate(est = purrr::map(q_daily_vec, estimate_q_hac)) %>%
    tidyr::unnest_wider(est) %>%
    mutate(
      q_lower = ifelse(is.finite(se_annual), q_hat_annual - 1.96 * se_annual, NA_real_),
      q_upper = ifelse(is.finite(se_annual), q_hat_annual + 1.96 * se_annual, NA_real_)
    ) %>%
    dplyr::select(State, N_days, q_hat_annual, se_q_annual = se_annual, q_lower, q_upper) %>%
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
    State       = states,
    sigma       = pars$sigma,
    se_sigma    = pars$se_sigma,
    sigma_lower = pars$sigma_lower,
    sigma_upper = pars$sigma_upper
  )

  mu_df %>%
    left_join(q_df,     by = "State") %>%
    left_join(sigma_df, by = "State") %>%
    mutate(
      Model  = model_label,
      mu_tot = mu_cap + q_hat_annual
    ) %>%
    arrange(State)
}

# ---- Formatting helpers (NA ⇒ dashes) ----
format_ci <- function(est, lower, upper, digits = 4) {
  fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)

  if (!is.finite(est)) {
    return("\\textemdash\\; (\\textemdash, \\textemdash)")
  }
  if (!is.finite(lower) || !is.finite(upper)) {
    return(paste0(fmt(est), "\\; (\\textemdash, \\textemdash)"))
  }
  paste0(fmt(est), "\\; (", fmt(lower), ", ", fmt(upper), ")")
}

format_se <- function(se, digits = 5) {
  if (!is.finite(se)) return("\\textemdash")
  sprintf(paste0("%.", digits, "f"), se)
}

build_param_rows <- function(state_df) {
  rows <- list()

  # mu_cap,i
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble(
      Parameter             = sprintf("$\\widehat{\\mu}_{\\text{cap},%d}$", s),
      `Estimate (95\\% CI)` = format_ci(r$mu_cap, r$mu_lower, r$mu_upper),
      `Std.\\ Error`        = format_se(r$se_mu)
    )
  }

  # q_i
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble(
      Parameter             = sprintf("$\\widehat{q}_{%d}$", s),
      `Estimate (95\\% CI)` = format_ci(r$q_hat_annual, r$q_lower, r$q_upper),
      `Std.\\ Error`        = format_se(r$se_q_annual)
    )
  }

  # mu_tot,i (no SE/CI)
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State

    est_str <- if (is.finite(r$mu_tot)) sprintf("%.4f", r$mu_tot) else "\\textemdash"

    rows[[length(rows) + 1]] <- tibble(
      Parameter             = sprintf("$\\widehat{\\mu}_{\\text{tot},%d}$", s),
      `Estimate (95\\% CI)` = paste0(est_str, "\\; (\\textemdash, \\textemdash)"),
      `Std.\\ Error`        = "\\textemdash"
    )
  }

  # sigma_i
  for (i in seq_len(nrow(state_df))) {
    r <- state_df[i, ]
    s <- r$State
    rows[[length(rows) + 1]] <- tibble(
      Parameter             = sprintf("$\\widehat{\\sigma}_{%d}$", s),
      `Estimate (95\\% CI)` = format_ci(r$sigma, r$sigma_lower, r$sigma_upper),
      `Std.\\ Error`        = format_se(r$se_sigma)
    )
  }

  bind_rows(rows)
}

# ---- LaTeX table builder ----
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
              param_rows$`Estimate (95\\% CI)`[i],
              param_rows$`Std.\\ Error`[i])
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

# ---- Model list (12 models) ----
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

# ---- Build tables for all models ----
latex_tables <- list()

for (k in seq_len(nrow(model_specs))) {
  spec <- model_specs[k, ]
  fit  <- models_BSHMM_2019[[ spec$name ]]
  stopifnot(!is.null(fit))

  message("Processing model: ", spec$name)

  state_df <- get_statewise_params_with_divs(
    fit         = fit,
    N           = spec$N,
    family      = spec$family,
    model_label = spec$name,
    spx_div     = spx_div,
    x_ret       = x_ret,
    dt          = dt
  )

  param_rows <- build_param_rows(state_df)

  cap <- sprintf(
    "%d-state BS-HMM with %s. State-dependent parameter estimates (annualised) with 95\\%% confidence intervals based on the inverse Hessian of the minimised log-likelihood where available. The capital-gains drift $\\mu_{\\text{cap},i}$ and dividend yield $q_i$ combine to give the total-return drift $\\mu_{\\text{tot},i} = \\mu_{\\text{cap},i} + q_i$ in each state. Standard errors and confidence intervals for $\\mu_{\\text{tot},i}$ are not reported and are indicated by dashes. Where the Hessian is singular or ill-conditioned, standard errors and confidence intervals for $\\mu_{\\text{cap},i}$ and $\\sigma_i$ are omitted and indicated by dashes. Where state-wise dividend observations are insufficient for a HAC estimate, standard errors and confidence intervals for $q_i$ are omitted and indicated by dashes.",
    spec$N,
    family_desc[[ spec$family ]]
  )

  lab <- sprintf("tab:BSHMM_%s", spec$name)

  latex_tables[[ spec$name ]] <- build_latex_table(
    param_rows = param_rows,
    caption    = cap,
    label      = lab
  )
}

# ---- Print all tables to console ----
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}

# ---- Write all tables to a single .tex file ----
out_file <- "BSHMM_state_tables_all_12_models_robust.tex"
sink(out_file)
for (nm in names(latex_tables)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables[[nm]], "\n\n")
}
sink()
cat("Wrote:", out_file, "\n")
