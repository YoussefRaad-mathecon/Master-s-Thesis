# BS Simulation


setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Parameters ----
n <- 25000
dt <- 1/252
mu_true <- 0.05
sigma_true <- 0.15
S0 <- 100

# ---- Simulation (exact GBM) ----
simulate_BS <- function(n, mu, sigma, S0, dt) {
  Z   <- rnorm(n)
  inc <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z
  logS <- c(log(S0), log(S0) + cumsum(inc))
  exp(logS)
}

# ---- Log-likelihood (log-returns) ----
BS_loglik <- function(params, S, dt) {
  mu <- params[1]
  sigma <- exp(params[2])
  R <- diff(log(S))
  mean_R <- (mu - 0.5 * sigma^2) * dt
  var_R  <- sigma^2 * dt
  ll <- sum(dnorm(R, mean = mean_R, sd = sqrt(var_R), log = TRUE))
  return(-ll)
}

# ---- Fit (nlm) ----
fit_BS <- function(S, dt) {
  init_params <- c(mu_true, log(sigma_true))
  result <- nlm(f = BS_loglik, p = init_params, S = S, dt = dt, print.level = 0)
  return(c(result$estimate[1], exp(result$estimate[2])))
}

# ---- Simulations ----
num_simulations <- 1

simulated_series_BS <- vector("list", num_simulations)
fitted_params_BS <- matrix(NA, nrow = num_simulations, ncol = 2)

for (i in 1:num_simulations) {
  set.seed(i)
  simulated_series_BS[[i]] <- simulate_BS(n, mu_true, sigma_true, S0, dt)
  fitted_params_BS[i, ] <- fit_BS(simulated_series_BS[[i]], dt)
}

BSFit <- numeric(2)
for (i in 1:2){
  BSFit[i] <- mean(fitted_params_BS[, i])
}

colnames(fitted_params_BS) <- c("mu", "sigma")
print(colMeans(fitted_params_BS))

save(fitted_params_BS, file = "fitted_params_BS.RData")

# ---- Plot one path ----
S <- simulated_series_BS[[1]]
t_years <- (0:(length(S)-1)) * dt

op <- options(scipen = 999)

library(ggplot2)

df <- data.frame(Years = t_years, StockPrice = S)

ggplot(df, aes(x = Years, y = StockPrice)) +
  geom_line(linewidth = 0.9, color = "#901a1E") +
  labs(x = "Years", y = "Stock Price") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

options(op)
