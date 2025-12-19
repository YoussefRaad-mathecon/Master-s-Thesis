# BS–SSM_beta: Hessian-based SEs for robustness grid + parameter tables

# ---- Paths ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
GRID_FILE <- "fitted_params_BSSSM_grid_2019_beta.RData"

# ---- Packages ----
library(dplyr)

# ---- Constant q from BSM ----
q_hat <- 0.03404597
q_se  <- 0.00092195
z     <- qnorm(0.975)
q_lower <- q_hat - z * q_se
q_upper <- q_hat + z * q_se

# ---- Load robustness grid (beta) ----
stopifnot(file.exists(GRID_FILE))
load(GRID_FILE)
stopifnot(exists("models_BSSSM_robustness_beta"))

grid_list_beta    <- models_BSSSM_robustness_beta$grid
grid_summary_beta <- models_BSSSM_robustness_beta$summary

# ---- Hessian inversion helper ----
invert_safe <- function(A) {
  if (any(!is.finite(A))) {
    message("Hessian contains non-finite entries; skipping SEs for this model.")
    return(matrix(NA_real_, nrow = nrow(A), ncol = ncol(A)))
  }
  
  tryCatch(
    solve(A),
    error = function(e) {
      message("Hessian not positive-definite; using pseudoinverse via MASS::ginv.")
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Please install.packages('MASS').")
      }
      MASS::ginv(A)
    }
  )
}

# ---- Natural-scale extraction (one grid entry) ----
extract_bsssm_beta_natural <- function(entry) {
  if (is.null(entry) || !is.list(entry) || entry$status != "ok") return(NULL)
  if (is.null(entry$fit)) return(NULL)
  
  fit <- entry$fit
  if (is.null(fit$estimate) || is.null(fit$hessian)) return(NULL)
  
  eta <- fit$estimate
  H_w <- fit$hessian
  
  if (!all(dim(H_w) == c(6, 6))) {
    stop("Unexpected Hessian dimension (expected 6x6) for m = ", entry$m, ", bm = ", entry$bm)
  }
  
  phi_hat        <- plogis(eta[1])
  sigZ_hat       <- exp(eta[2])
  mu0_hat        <- eta[3]
  sigma_hat      <- exp(eta[4])
  beta_mu_hat    <- eta[5]
  beta_sigma_hat <- eta[6]
  
  theta_hat <- c(
    phi        = phi_hat,
    sigZ       = sigZ_hat,
    mu_cap     = mu0_hat,
    sigma      = sigma_hat,
    beta_mu    = beta_mu_hat,
    beta_sigma = beta_sigma_hat
  )
  
  M <- diag(c(
    phi_hat * (1 - phi_hat),
    sigZ_hat,
    1,
    sigma_hat,
    1,
    1
  ))
  rownames(M) <- colnames(M) <- names(theta_hat)
  
  Cov_w <- invert_safe(H_w)
  Cov_n <- M %*% Cov_w %*% t(M)
  
  var_vec <- pmax(diag(Cov_n), 0)
  se      <- sqrt(var_vec)
  
  z <- 1.96
  CI_lower <- theta_hat - z * se
  CI_upper <- theta_hat + z * se
  
  list(
    theta_hat = theta_hat,
    se        = se,
    CI_lower  = CI_lower,
    CI_upper  = CI_upper
  )
}

# ---- LaTeX formatting helpers ----
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

format_se <- function(se, digits = 6) {
  if (!is.finite(se)) return("\\textemdash")
  sprintf(paste0("%.", digits, "f"), se)
}

build_bsssm_beta_rows <- function(theta_hat, se, CI_lower, CI_upper,
                                  q_hat, q_se, q_lower, q_upper) {
  rows <- list()
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\rho}$",
    Est_CI    = format_ci(theta_hat["phi"], CI_lower["phi"], CI_upper["phi"]),
    SE        = format_se(se["phi"])
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\sigma}_{\\varepsilon}$",
    Est_CI    = format_ci(theta_hat["sigZ"], CI_lower["sigZ"], CI_upper["sigZ"]),
    SE        = format_se(se["sigZ"])
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{q}$",
    Est_CI    = format_ci(q_hat, q_lower, q_upper),
    SE        = format_se(q_se)
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\mu}_{\\text{cap}}$",
    Est_CI    = format_ci(theta_hat["mu_cap"], CI_lower["mu_cap"], CI_upper["mu_cap"]),
    SE        = format_se(se["mu_cap"])
  )
  
  mu_tot_hat <- theta_hat["mu_cap"] + q_hat
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\mu}_{\\text{tot}}$",
    Est_CI    = paste0(sprintf("%.4f", mu_tot_hat), "\\; (\\textemdash, \\textemdash)"),
    SE        = "\\textemdash"
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\sigma}$",
    Est_CI    = format_ci(theta_hat["sigma"], CI_lower["sigma"], CI_upper["sigma"]),
    SE        = format_se(se["sigma"])
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\beta}_{\\mu}$",
    Est_CI    = format_ci(theta_hat["beta_mu"], CI_lower["beta_mu"], CI_upper["beta_mu"]),
    SE        = format_se(se["beta_mu"])
  )
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    Parameter = "$\\widehat{\\beta}_{\\sigma}$",
    Est_CI    = format_ci(theta_hat["beta_sigma"], CI_lower["beta_sigma"], CI_upper["beta_sigma"]),
    SE        = format_se(se["beta_sigma"])
  )
  
  dplyr::bind_rows(rows)
}

build_latex_table_bsssm_beta <- function(param_rows, caption, label) {
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\begin{tabular}{lrr}",
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
    "\\captionsetup{type=table}",
    sprintf("\\captionof{table}{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\end{table}"
  )
  
  paste(lines, collapse = "\n")
}

# ---- Build LaTeX tables across grid (beta) ----
latex_tables_beta <- list()

for (i in seq_along(grid_list_beta)) {
  entry <- grid_list_beta[[i]]
  
  if (is.null(entry) || !is.list(entry) || entry$status != "ok") {
    message("Skipping combo ", i, " (status != 'ok').")
    next
  }
  
  nat <- extract_bsssm_beta_natural(entry)
  if (is.null(nat)) {
    message("Skipping combo ", i, " (no natural-scale params).")
    next
  }
  
  theta_hat <- nat$theta_hat
  se        <- nat$se
  CI_lower  <- nat$CI_lower
  CI_upper  <- nat$CI_upper
  
  param_rows <- build_bsssm_beta_rows(
    theta_hat = theta_hat,
    se        = se,
    CI_lower  = CI_lower,
    CI_upper  = CI_upper,
    q_hat     = q_hat,
    q_se      = q_se,
    q_lower   = q_lower,
    q_upper   = q_upper
  )
  
  m_val  <- entry$m
  bm_val <- entry$bm
  name   <- names(grid_list_beta)[i]
  
  cap <- sprintf(
    "BS-SSM$_{\\beta}$ using an $m = %d$ point grid and truncation $b_{\\max} = %.1f$. Parameter estimates with 95\\%% confidence intervals in parentheses.",
    m_val, bm_val
  )
  
  lab <- sprintf("tab:BSSSMbeta_m%d_bm%.1f", m_val, bm_val)
  lab <- gsub("\\.", "", lab)
  
  latex_tables_beta[[name]] <- build_latex_table_bsssm_beta(
    param_rows,
    caption = cap,
    label   = lab
  )
}

# ---- Print all tables ----
for (nm in names(latex_tables_beta)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables_beta[[nm]], "\n\n")
}

# ---- Write all tables to file ----
sink("BSSSM_param_tables_grid_2019_beta.tex")
for (nm in names(latex_tables_beta)) {
  cat("% ------------------------------------------------------------------\n")
  cat("%", nm, "\n")
  cat(latex_tables_beta[[nm]], "\n\n")
}
sink()
