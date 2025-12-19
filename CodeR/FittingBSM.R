####################################################################################################################
####################################################################################################################
#---------------------------------------- Fitting BS ----------------------------------------------------------
####################################################################################################################
##################################################################################################################

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
set.seed(112233)
load("sp500_clean.Rdata")


### Parameters
dt <- 1/252
spx <- spx %>% slice_head(n = 23109)     

S <- na.omit(spx$Return)
#----------------------------------------------Log-likelihood--------------------------------------------------#
BS_loglik <- function(params, S) {
  mu <- params[1]
  sigma <- exp(params[2])  # keep positivity via log-param like your Vasicek code
  
  n <- length(S)
  # log returns
  R <- na.omit(spx$Return)
  
  mean_R <- (mu - 0.5 * sigma^2) * dt
  var_R  <- sigma^2 * dt
  
  ll <- sum(dnorm(R, mean = mean_R, sd = sqrt(var_R), log = TRUE))
  return(-ll)  # negative log-likelihood for minimization
}

#-------------------------------------------Fitting---------------------------------------------------------------------#
# Continue with the BS log-likelihood optimization as before
mod <- function(S) {
  init_params <- c(0.06, log(0.2))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)
#------------------------------------------------------------------------------------------------------#
mod <- function(S) {
  init_params <- c(0.02, log(0.6))
  result <- nlm(f = BS_loglik, p = init_params, S = S, print.level = 2)
  return(c(result$estimate[1], exp(result$estimate[2])))
}
mod(S)
#------------------------------------------------------------------------------------------------------#
# one-line fit (matches your nlm HMM call style)
mod_BS <- nlm(BS_loglik, c(0.06, log(0.2)), S = S, print.level = 2, iterlim = 10000, hessian = TRUE)

# one-line print (k = 2)
cat(sprintf("mu=%.6f  sigma=%.6f  AIC=%.3f  BIC=%.3f\n",
            mod_BS$estimate[1], exp(mod_BS$estimate[2]),
            2*2 + 2*mod_BS$minimum,
            2*log(length(S)) + 2*mod_BS$minimum))


models_BS_2019 <- list(fit_BS)

save(models_BS_2019, file = "fitted_params_BS_2019.RData")


fit_BS
# Pieces for AIC/BIC
k          <- 2                         # parameters: mu and sigma
n_eff      <- length(S)                 # sample size used in the likelihood
logLik_hat <- -fit_BS$minimum              # because BS_loglik returns NEGATIVE log-lik

# Compute + print
AIC  <- 2*k - 2*logLik_hat
BIC  <- k*log(n_eff) - 2*logLik_hat

cat(sprintf("AIC: %.3f\nBIC: %.3f\n", AIC, BIC))


################################################################################################
#--------------------------- Closed-form MLEs for Black–Scholes (no optimization) --------------#
################################################################################################
bs_mle <- function(R, dt) {
  R   <- as.numeric(na.omit(R))
  Tn  <- length(R)
  xbar <- mean(R)
  
  # MLEs from your derivation:
  sigma2_hat <- sum((R - xbar)^2) / (Tn * dt)     # \hat{\sigma}^2
  mu_hat     <- xbar / dt + 0.5 * sigma2_hat      # \hat{\mu}
  sigma_hat  <- sqrt(sigma2_hat)
  
  # Log-likelihood at the MLE:
  var_hat    <- sigma2_hat * dt                   # = \hat{\sigma}^2 \Delta
  logLik_hat <- -(Tn/2) * (log(2*pi*var_hat) + 1) # = -T/2 [ log(2π σ^2 Δ) + 1 ]
  
  # Information criteria (k = 2: mu and sigma)
  k   <- 2L
  AIC <- 2*k - 2*logLik_hat
  BIC <- k*log(Tn) - 2*logLik_hat
  
  list(mu = mu_hat,
       sigma = sigma_hat,
       sigma2 = sigma2_hat,
       T = Tn,
       dt = dt,
       xbar = xbar,
       logLik = logLik_hat,
       AIC = AIC,
       BIC = BIC)
}

# Run closed-form MLEs on your data:
fit_BS_mle <- bs_mle(S, dt)

cat(sprintf("Closed-form MLE  ->  mu=%.8f  sigma=%.8f  logLik=%.3f  AIC=%.3f  BIC=%.3f\n",
            fit_BS_mle$mu, fit_BS_mle$sigma, fit_BS_mle$logLik, fit_BS_mle$AIC, fit_BS_mle$BIC))

