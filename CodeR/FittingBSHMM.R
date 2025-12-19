# BS–HMM (S&P 500 to 2019): Fit all model variants (μ, σ, μσ families, 2–5 states)

# ---- Packages ----
library(tidyverse)

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- Data ----
load("sp500_clean.Rdata")
spx <- spx %>% slice_head(n = 23109)     
S <- na.omit(spx$Return)
x_raw <- S           # already log returns
x     <- as.numeric(x_raw)   # coerce to numeric vector


### Parameters
dt      <- 1/252


#------------------------------------------ mu ------------------------------------------------------------


### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N, dt = 1/252){
  mu <- theta.star[1:N]
  sigma <- exp(theta.star[N+1])
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(N+2):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  #delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  Tn <- length(x)
  allprobs <- matrix(1, Tn, N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[, j] <- dnorm(x[ind], mean = (mu[j] - 0.5 * sigma[1]^2) * dt, sd = sigma[1] * sqrt(dt), log = FALSE)
    
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



### 2-state HMM
N <- 2
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(0.1),
                rep(-2, (N - 1) * N))
mod2_mu <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 3-state HMM
N <- 3
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(0.1),
                rep(-2, (N - 1) * N))
mod3_mu <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)

### 4-state HMM
N <- 4
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(0.1),
                rep(-2, (N - 1) * N))
mod4_mu <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 5-state HMM
N <- 5
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(0.1),
                rep(-2, (N - 1) * N))
mod5_mu <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)





#------------------------------------------ sigma ------------------------------------------------------------

### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N, dt = 1/252){
  mu <- theta.star[1]
  sigma <- exp(theta.star[2:(N+1)])
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(N+2):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  #delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  Tn <- length(x)
  allprobs <- matrix(1, Tn, N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[, j] <- dnorm(x[ind], mean = (mu[1] - 0.5 * sigma[j]^2) * dt, sd = sigma[j] * sqrt(dt), log = FALSE)
    
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


### 2-state HMM
N <- 2
theta.star <- c((0.02),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod2_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)



### 3-state HMM
N <- 3
theta.star <- c((0.02),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod3_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 4-state HMM
N <- 4
theta.star <- c((0.02),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod4_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 5-state HMM
N <- 5
theta.star <- c((0.02),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod5_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)



#------------------------------------------ mu & sigma ------------------------------------------------------------


### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N, dt = 1/252){
  mu <- theta.star[1:N]
  sigma <- exp(theta.star[(N+1):(N*2)])
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[((N*2)+1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  #delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  Tn <- length(x)
  allprobs <- matrix(1, Tn, N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[, j] <- dnorm(x[ind], mean = (mu[j] - 0.5 * sigma[j]^2) * dt, sd = sigma[j] * sqrt(dt), log = FALSE)
    
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


### 2-state HMM
N <- 2
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod2_mu_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)



### 3-state HMM
N <- 3
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod3_mu_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 4-state HMM
N <- 4
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod4_mu_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)


### 5-state HMM
N <- 5
theta.star <- c((seq(0.02, 0.06, length = N)),
                log(seq(0.1, 0.2, length = N)),
                rep(-2, (N - 1) * N))
mod5_mu_sigma <- nlm(mllk, theta.star, x = x, N = N, dt = 1/252, print.level = 2, iterlim = 10000, hessian = TRUE)

### Extracting tpm
Gamma5_hat <- diag(N)
Gamma5_hat[!Gamma5_hat] <- exp(mod5_mu_sigma$estimate[((N*2)+1):(length(theta.star))])
Gamma5_hat <- Gamma5_hat / rowSums(Gamma5_hat)
Gamma5_hat

### Extracting stationary distribution
delta5_hat <- solve(t(diag(N) - Gamma5_hat + 1), rep(1, N))
delta5_hat


## After running nlm(...)
n  <- length(x)

## Parameter count used in your parameterization:
## N means + N log(sigmas) + N*(N-1) transition logits  =>  k = N^2 + N
k  <- N^2 + N

## nlm minimized mllk, which is (usually) -logLik:
logLik_hat <- -mod7_mu_sigma$minimum   # <-- if mllk is +logLik, use: logLik_hat <- mod7_mu_sigma$minimum

AIC_val <- 2*k - 2*logLik_hat
BIC_val <- log(n)*k - 2*logLik_hat

c(logLik = logLik_hat, k = k, AIC = AIC_val, BIC = BIC_val)




models_BSHMM_2019 <- list(
  mod2_mu = mod2_mu,
  mod3_mu = mod3_mu,
  mod4_mu = mod4_mu,
  mod5_mu = mod5_mu,
  
  mod2_sigma = mod2_sigma,
  mod3_sigma = mod3_sigma,
  mod4_sigma = mod4_sigma,
  mod5_sigma = mod5_sigma,
  
  mod2_mu_sigma = mod2_mu_sigma,
  mod3_mu_sigma = mod3_mu_sigma,
  mod4_mu_sigma = mod4_mu_sigma,
  mod5_mu_sigma = mod5_mu_sigma
)

save(models_BSHMM_2019, file = "fitted_params_BSHMM_2019.RData")
load("fitted_params_BSHMM_2019.RData")
