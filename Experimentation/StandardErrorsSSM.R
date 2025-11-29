################################################################################
# Black–Scholes SSM (GBM): Load fit and compute Hessian-based standard errors
################################################################################

## --- User settings ------------------------------------------------------------
# Set your working directory to where both files live:
# setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

DATA_FILE <- "sp500_clean.Rdata"                # contains 'spx' (your cleaned SP500 df)
FIT_FILE  <- "fitted_params_BSSSM_2019.RData"   # contains list 'models_BSSSM_2019'
OUT_CSV   <- "BSSSM_SEs.csv"                    # optional output

## --- Load data frame (optional, for inspection) -------------------------------
if (file.exists(DATA_FILE)) {
  load(DATA_FILE)   # should create 'spx'
  if (exists("spx")) {
    cat("Data frame 'spx' loaded.\n")
    cat("Columns:", paste(colnames(spx), collapse = ", "), "\n")
    cat("Rows   :", nrow(spx), "\n\n")
  } else {
    warning("Loaded '", DATA_FILE, "' but did not find object 'spx'.")
  }
} else {
  warning("Data file '", DATA_FILE, "' not found. Skipping data load.")
}

## --- Load saved fit -----------------------------------------------------------
stopifnot(file.exists(FIT_FILE))
load(FIT_FILE)  # expects 'models_BSSSM_2019'
stopifnot(exists("models_BSSSM_2019"), length(models_BSSSM_2019) >= 1)

# If you have multiple fits saved, pick one here (default: first)
fit_BS <- models_BSSSM_2019[[1]]

# Sanity checks
stopifnot(is.list(fit_BS))
stopifnot(!is.null(fit_BS$estimate), !is.null(fit_BS$hessian))
H_w <- fit_BS$hessian
theta_hat_tr <- fit_BS$estimate  # working-scale MLEs from 'nlm'

if (!all(dim(H_w) == c(4,4))) stop("Unexpected Hessian dimension (expected 4x4).")

## --- Back-transform to natural scale ------------------------------------------
# Working params used when fitting:
#   eta1 = qlogis(phi), eta2 = log(sigZ), eta3 = mu, eta4 = log(sigma)
theta_hat <- c(
  phi   = plogis(theta_hat_tr[1]),
  sigZ  = exp(theta_hat_tr[2]),
  mu    = theta_hat_tr[3],
  sigma = exp(theta_hat_tr[4])
)

## --- Build Jacobian M = d(theta)/d(eta) at the MLE ----------------------------
# d phi / d qlogis(phi) = phi*(1-phi)
# d sigZ / d log(sigZ)  = sigZ
# d mu  / d mu          = 1
# d sigma / d log(sigma)= sigma
M <- diag(c(
  theta_hat["phi"]  * (1 - theta_hat["phi"]),
  theta_hat["sigZ"],
  1,
  theta_hat["sigma"]
))
rownames(M) <- colnames(M) <- names(theta_hat)

## --- Safe inversion helper ----------------------------------------------------
invert_safe <- function(A) {
  tryCatch(solve(A), error = function(e) {
    message("Hessian not positive-definite; using pseudoinverse via MASS::ginv.")
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install.packages('MASS').")
    MASS::ginv(A)
  })
}

## --- Covariance on working scale, then natural scale --------------------------
# nlm returns Hessian of the objective (-log L). Under regularity,
# Var_hat(eta_hat) ≈ H_w^{-1}. With H_w = M H_n M^T, we get H_n^{-1} = M H_w^{-1} M^T.
Cov_w <- invert_safe(H_w)
Cov_n <- M %*% Cov_w %*% t(M)

## --- Standard errors and Wald CIs on natural scale ----------------------------
se <- sqrt(pmax(diag(Cov_n), 0))  # guard against tiny negative due to numerics
z  <- 1.96
CI_lower <- theta_hat - z * se
CI_upper <- theta_hat + z * se

se_table <- data.frame(
  Parameter = names(theta_hat),
  Estimate  = as.numeric(theta_hat),
  StdError  = as.numeric(se),
  CI95_L    = as.numeric(CI_lower),
  CI95_U    = as.numeric(CI_upper),
  row.names = NULL,
  check.names = FALSE
)

## --- print
cat("\nParameter estimates with Hessian-based SEs (natural scale):\n")
print(transform(se_table,
                Estimate = signif(Estimate, 6),
                StdError = signif(StdError, 6),
                CI95_L   = signif(CI95_L, 6),
                CI95_U   = signif(CI95_U, 6)
), row.names = FALSE)

