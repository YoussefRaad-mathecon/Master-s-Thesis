####################################################################################################################
####################################################################################################################
#---------------------------------------- BS Simulation ------------------------------------------------------------
####################################################################################################################
####################################################################################################################
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

### Set parameters for the Black–Scholes process
n <- 25000
dt <- 1/252
mu_true <- 0.05
sigma_true <- 0.15
S0 <- 100

### Function to simulate Black–Scholes (GBM) prices (exact discretization)
simulate_BS <- function(n, mu, sigma, S0, dt) {
  # returns length n+1, including S0
  Z   <- rnorm(n)
  inc <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z
  logS <- c(log(S0), log(S0) + cumsum(inc))
  exp(logS)
}

### Log-likelihood for GBM using log returns
# If R_i = log(S_i) - log(S_{i-1}), then R_i ~ N( (mu - 0.5*sigma^2)*dt, sigma^2 * dt )
BS_loglik <- function(params, S, dt) {
  mu <- params[1]
  sigma <- exp(params[2])  # positivity via log-param
  R <- diff(log(S))        # *** FIX: use log-returns ***
  mean_R <- (mu - 0.5 * sigma^2) * dt
  var_R  <- sigma^2 * dt
  ll <- sum(dnorm(R, mean = mean_R, sd = sqrt(var_R), log = TRUE))
  return(-ll)  # negative log-likelihood for minimization
}

### Function to fit the BS model (matches your nlm pattern & transforms)
fit_BS <- function(S, dt) {
  init_params <- c(mu_true, log(sigma_true))
  result <- nlm(f = BS_loglik, p = init_params, S = S, dt = dt, print.level = 0)
  return(c(result$estimate[1], exp(result$estimate[2])))
}

### Number of simulations
num_simulations <- 1

### Storage
simulated_series_BS <- vector("list", num_simulations)
fitted_params_BS <- matrix(NA, nrow = num_simulations, ncol = 2)

### Simulate and fit
for (i in 1:num_simulations) {
  set.seed(i)  # reproducibility
  
  # *** FIX: pass args and use multiplicative simulator ***
  simulated_series_BS[[i]] <- simulate_BS(n, mu_true, sigma_true, S0, dt)
  
  # Fit BS on the simulated path
  fitted_params_BS[i, ] <- fit_BS(simulated_series_BS[[i]], dt)
}

BSFit <- numeric(2)
for (i in 1:2){
  BSFit[i] <- mean(fitted_params_BS[, i])
}

# Quick summary
colnames(fitted_params_BS) <- c("mu", "sigma")
print(colMeans(fitted_params_BS))

save(fitted_params_BS, file = "fitted_params_BS.RData")

# Pick one simulated path (first)
S <- simulated_series_BS[[1]]

# Time in years (matches length of S)
t_years <- (0:(length(S)-1)) * dt

# ggplot
op <- options(scipen = 999)  # avoid scientific notation just in case
library(ggplot2)

df <- data.frame(Years = t_years, StockPrice = S)

ggplot(df, aes(x = Years, y = StockPrice)) +
  geom_line(linewidth = 0.9, color = "#901a1E") +
  labs(x = "Years", y = "Stock Price") +
  theme_bw() +  # keeps gridlines
  theme(
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

options(op)
