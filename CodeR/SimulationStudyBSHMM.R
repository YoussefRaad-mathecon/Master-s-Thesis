# BS–HMM: Simulation + MLE on log-returns


# ---- Path ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Packages ----
library(ggplot2)
library(gridExtra)
library(GMCM)
library(matrixStats)
library(scales)


# ---- Parameters ----
n  <- 25000
dt <- 1/252
N  <- 5

mu_true    <- c(-0.20, -0.06, 0.03, 0.08, 0.12)
sigma_true <- c( 0.30,  0.24, 0.18, 0.15, 0.22)
S0 <- 100

Gamma_true <- diag(N)
Gamma_true[!Gamma_true] <- seq(0.01, 0.05, length = N * (N - 1))
Gamma_true <- Gamma_true / rowSums(Gamma_true)
delta_true <- solve(t(diag(N) - Gamma_true + 1), rep(1, N))

Gamma_true
delta_true

# ---- Simulation (regime-switching GBM; returns log-returns) ----
simulate_HMM_BS <- function(n, N, mu_true, sigma_true, Gamma_true, delta_true, S0, dt = 1/252) {
  State <- numeric(n)
  S     <- numeric(n)
  State[1] <- sample(x = 1:N, size = 1, prob = delta_true)
  S[1] <- log(S0)
  
  for (i in 2:n) {
    State[i] <- sample(x = 1:N, size = 1, prob = Gamma_true[State[i - 1], ])
    
    j <- State[i]
    Z <- rnorm(1)
    S[i] <- (S[i-1]) + (mu_true[j] - 0.5 * sigma_true[j]^2) * dt + sigma_true[j] * Z * sqrt(dt)
  }
  
  r <- diff(S)
  return(list(r = r, State = State, S = S))
}

# ---- Log-likelihood for BS–HMM (log-returns) ----
mllk_BS <- function(theta.star, x, N, dt = 1/252) {
  mu    <- theta.star[1:N]
  sigma <- exp(theta.star[(N + 1):(2 * N)])
  
  Gamma <- diag(N)
  off_diag <- exp(theta.star[(2 * N + 1):length(theta.star)])
  Gamma[!Gamma] <- off_diag
  Gamma <- Gamma / rowSums(Gamma)
  
  delta <- rep(1 / N, N)
  
  Tn <- length(x)
  allprobs <- matrix(1, Tn, N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[, j] <- dnorm(x[ind],
                           mean = (mu[j] - 0.5 * sigma[j]^2) * dt,
                           sd   = sigma[j] * sqrt(dt),
                           log  = FALSE)
  }
  
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)){
    if (t > nrow(allprobs)) {
      stop(paste("Subscript out of bounds: t =", t, "exceeds nrow(allprobs) =", nrow(allprobs)))
    }
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  
  return(-l)
}

# ---- Fit (nlm) ----
fit_HMM_BS <- function(r, N, theta.star, dt = 1/252) {
  result <- nlm(mllk_BS, theta.star, x = r, N = N, dt = dt, print.level = 0, iterlim = 10000)
  return(c(
    result$estimate[1:N],
    exp(result$estimate[(N + 1):(2 * N)]),
    exp(result$estimate[(2 * N + 1):length(theta.star)])
  ))
}

# ---- Simulate and fit ----
num_simulations <- 1

theta.star <- c(
  mu_true,
  log(sigma_true),
  log(Gamma_true[row(Gamma_true) != col(Gamma_true)])
)

fitted_params_HMM <- matrix(NA, nrow = num_simulations, ncol = length(theta.star))

set.seed(1)
for (i in 1:num_simulations) {
  cat(i, "\n")
  sim <- simulate_HMM_BS(n, N, mu_true, sigma_true, Gamma_true, delta_true, S0, dt = dt)
  r <- sim$r
  fitted_params_HMM[i, ] <- fit_HMM_BS(r, N, theta.star, dt = dt)
}

print(fitted_params_HMM)

# ---- Means across simulations ----
HMMFit <- numeric(length(theta.star))
for (i in 1:length(theta.star)) HMMFit[i] <- mean(fitted_params_HMM[, i])

Gamma5Sim <- diag(N)
Gamma5Sim[!Gamma5Sim] <- HMMFit[(2 * N + 1):length(theta.star)]
Gamma5Sim <- Gamma5Sim / rowSums(Gamma5Sim)

solve(t(diag(N) - Gamma5Sim + 1), rep(1, N))

# ---- Summaries / scaffolding ----
SimHMM <- as.data.frame(fitted_params_HMM)
colnames(SimHMM) <- c(
  paste0("mu", 1:N),
  paste0("sigma", 1:N),
  "P(2->1)", "P(3->1)", "P(4->1)", "P(5->1)",
  "P(1->2)", "P(3->2)", "P(4->2)", "P(5->2)",
  "P(1->3)", "P(2->3)", "P(4->3)", "P(5->3)",
  "P(1->4)", "P(2->4)", "P(3->4)", "P(5->4)",
  "P(1->5)", "P(2->5)", "P(3->5)", "P(4->5)"
)
SimHMM$Index <- 1:num_simulations
colMeans(SimHMM[, 1:(ncol(SimHMM)-1)])
apply(SimHMM[, 1:(ncol(SimHMM)-1)], 2, sd)

save(fitted_params_HMM, file = "fitted_params_BSHMM.RData")

# ---- Plot one simulated price path ----
S <- exp(sim$S)
t_years <- (0:(length(S)-1)) * dt

op <- options(scipen = 999)

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
