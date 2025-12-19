# BS–SSM vs BS–SSM_beta: AIC/BIC comparison for best grid models

# ---- Packages ----
library(dplyr)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- User settings ----
FILE_BSSSM      <- "BSSSM_grid_summary_sorted_2019.rds"
FILE_BSSSM_beta <- "BSSSM_grid_summary_sorted_2019_beta.rds"

# Best (m, bm) combinations (from your robustness tables)
best_BSSSM_m   <- 70
best_BSSSM_bm  <- 3
best_BSSSMb_m  <- 100
best_BSSSMb_bm <- 4

# Number of parameters in each model
p_BSSSM      <- 4  # rho, sigma_eps, mu_cap, sigma
p_BSSSM_beta <- 6  # rho, sigma_eps, mu0, sigma, beta_mu, beta_sigma

# ---- Load data & sample size (for BIC) ----
load("sp500_clean.Rdata")  # creates 'spx'
spx <- spx %>% slice_head(n = 23109)
x   <- as.numeric(na.omit(spx$Return))
n   <- length(x)

# ---- Load robustness summaries ----
stopifnot(file.exists(FILE_BSSSM), file.exists(FILE_BSSSM_beta))
grid_bs      <- readRDS(FILE_BSSSM)
grid_bs_beta <- readRDS(FILE_BSSSM_beta)

# ---- Extract negative log-likelihoods for chosen (m, bm) ----
get_negloglik <- function(df, m_val, bm_val) {
  row <- df %>% dplyr::filter(m == m_val, bm == bm_val)
  
  if (nrow(row) == 0L) stop("No row found for m = ", m_val, ", bm = ", bm_val)
  if (!"negloglik" %in% names(row)) stop("Column 'negloglik' not found in the data frame.")
  
  as.numeric(row$negloglik[1])
}

negloglik_BSSSM      <- get_negloglik(grid_bs,      best_BSSSM_m,  best_BSSSM_bm)
negloglik_BSSSM_beta <- get_negloglik(grid_bs_beta, best_BSSSMb_m, best_BSSSMb_bm)

# ---- Compute AIC and BIC ----
AIC_BSSSM      <- 2 * p_BSSSM      + 2 * negloglik_BSSSM
BIC_BSSSM      <- log(n) * p_BSSSM      + 2 * negloglik_BSSSM

AIC_BSSSM_beta <- 2 * p_BSSSM_beta + 2 * negloglik_BSSSM_beta
BIC_BSSSM_beta <- log(n) * p_BSSSM_beta + 2 * negloglik_BSSSM_beta

# ---- Build comparison table ----
fmt1 <- function(z) sprintf("%.1f", z)
bold <- function(s) paste0("\\textbf{", s, "}")

tbl <- tibble::tibble(
  Model = c("BS-SSM", "BS-SSM$_{\\beta}$"),
  p     = c(p_BSSSM, p_BSSSM_beta),
  AIC   = c(AIC_BSSSM, AIC_BSSSM_beta),
  BIC   = c(BIC_BSSSM, BIC_BSSSM_beta)
) %>%
  dplyr::mutate(
    AIC = ifelse(AIC == min(AIC), bold(fmt1(AIC)), fmt1(AIC)),
    BIC = ifelse(BIC == min(BIC), bold(fmt1(BIC)), fmt1(BIC))
  )

# ---- Output ----
print(tbl, n = Inf)

# ---- LaTeX table output ----
# library(knitr); library(kableExtra)
# tbl %>%
#   knitr::kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "AIC and BIC for the BS-SSM and BS-SSM$_{\\beta}$ fitted to daily S\\&P 500 log-returns. Here, $p$ denotes the number of estimated parameters. The globally lowest (best) AIC and BIC are shown in bold.",
#     label    = "tab:SSM2019AICBIC"
#   ) %>%
#   kableExtra::kable_styling(position = "center", latex_options = "hold_position")
