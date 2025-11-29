################################################################################################
#--------------------------- Fitting SSM: Black–Scholes (GBM) ---------------------------------
################################################################################################
### Load library, data and WD
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library("fHMM") ### Model fitting
library(Quandl) ### Data
library(dplyr) ### Data manipulations
library(lubridate) ### Time-variables
library(stats) ### ACF plots
library(matrixcalc) ### Matrix calculations
library("RColorBrewer") ### Colors
library(latex2exp) ### Text for plots
library(matrixStats) ### ColSds
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)
load("sp500_clean.Rdata")
spx <- spx %>% slice_head(n = 23109)   
S <- na.omit(spx$Return)
x_raw <- S           # already log returns
x     <- as.numeric(x_raw)   # coerce to numeric vector


# If not already in memory:
if (!exists("models_BSHMM_2019")) load("fitted_params_BSHMM_2019.RData")

# n = #observations used in the fit
n_obs <- if (exists("x")) length(x) else if (exists("S")) length(S) else stop("Set n_obs = length(returns).")

library(tidyverse)
library(stringr)

summarise_mod <- function(name, mod, n_obs){
  # parse N and family from name like "mod5_mu_sigma", "mod3_sigma", "mod2_mu"
  m <- str_match(name, "^mod(\\d+)_(mu_sigma|mu|sigma)$")
  N <- as.integer(m[2]); family <- m[3]
  
  k <- length(mod$estimate)         # total free parameters actually estimated
  negloglik <- mod$minimum          # -log L at optimum
  aic <- 2*k + 2*negloglik
  bic <- k*log(n_obs) + 2*negloglik
  
  # expected split (for a quick check)
  k_emit <- switch(family,
                   "mu"       = N + 1L,  # N means + 1 shared sigma
                   "sigma"    = 1L + N,  # 1 shared mean + N sigmas
                   "mu_sigma" = 2L*N)    # N means + N sigmas
  k_gamma <- N*(N-1L)               # Γ has N(N-1) free params (row sums = 1)
  tibble(
    model = name, family, N,
    negloglik, k, k_emit, k_gamma, k_expected = k_emit + k_gamma,
    k_ok = (k == k_expected),
    aic, bic,
    conv_code = mod$code
  )
}

summary_tbl <- purrr::imap_dfr(models_BSHMM_2019, ~ summarise_mod(.y, .x, n_obs)) %>%
  group_by(family) %>%
  mutate(
    aic_rank = rank(aic, ties.method = "min"),
    bic_rank = rank(bic, ties.method = "min"),
    delta_aic = aic - min(aic),
    delta_bic = bic - min(bic),
    aic_weight = exp(-0.5*delta_aic) / sum(exp(-0.5*delta_aic))
  ) %>%
  ungroup() %>%
  arrange(family, N)

# Full table
print(summary_tbl, n = Inf)

# Compact view (sorted by BIC within family)
summary_tbl %>%
  arrange(family, bic) %>%
  select(family, N, k, negloglik, aic, bic, aic_rank, bic_rank, aic_weight, k_ok, conv_code) %>%
  print(n = Inf)
# Exports if you like:
# write.csv(summary_tbl, "BSHMM_AICBIC_2019.csv", row.names = FALSE)
# save(summary_tbl, file = "BSHMM_AICBIC_2019.RData")

# (Optional LaTeX)
# library(knitr); library(kableExtra)
# summary_tbl %>%
#   arrange(family, bic) %>%
#   mutate(across(c(negloglik,aic,bic), ~round(.x, 2))) %>%
#   select(Family = family, States = N, k, `-log L` = negloglik, AIC = aic, BIC = bic,
#          `ΔAIC` = delta_aic, `Weight` = aic_weight) %>%
#   kable("latex", booktabs = TRUE, caption = "AIC/BIC for BS–HMM (S&P500 to 2019)") %>%
#   kable_styling(position = "center", latex_options = "hold_position")
