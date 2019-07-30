# --- simulate some data
set.seed(10)

t <- seq(1:300)
N <- length(t)
X <- 60 + truncnorm::rtruncnorm(N, 0, Inf, 10, 40)
e <- rep(4, N)

plot(t, X, ylim = c(0, max(X)))

#  -- the parameters
#   - beta = mean travel time
#   - kappa = between-vehicle variability
#   - epsilon, phi = system noise of beta, kappa respectively

beta <- array(NA, c(1, N + 1))
beta_hat <- beta
P <- array(NA, c(1, 1, N + 1))
P_hat <- P

kappa <- kappa_hat <- array(NA, c(1, N + 1))
R <- R_hat <- array(NA, c(1, 1, N + 1))

# initialize
beta[,1] <- 60
P[,,1] <- 100
kappa[,1] <- 30
R[,,1] <- 100

epsilon <- 0.2
phi <- 0.5
r <- 5

for (i in 1:N) {
    beta_hat[,i+1] <- beta[,i]
    P_hat[,,i+1] <- P[,,i] + epsilon^2

    kappa_hat[,i+1] <- kappa[,i]
    R_hat[,,i+1] <- R[,,i] + phi^2

    # first, update kappa
    z <- abs(X[i] - beta_hat[,i+1]) / e[i] - kappa_hat[,i+1, drop = FALSE]
    S <- as.matrix(R_hat[,,i+1, drop = FALSE] + r^2)
    K <- R_hat[,,i+1,drop=F] %*% solve(S)
    kappa[,i+1] <- kappa_hat[,i+1] + K %*% z
    R[,,i+1] <- (1 - K) %*% R_hat[,,i+1]

    y <- X[i] - beta_hat[,i+1, drop = FALSE]
    S <- as.matrix(P_hat[,,i+1, drop = FALSE] + e[i]^2)
    K <- P_hat[,,i+1, drop = FALSE] %*% solve(S)
    beta[,i+1] <- beta_hat[,i+1] + K %*% y
    P[,,i+1] <- (1 - K) %*% P_hat[,,i+1]
}

plot(t, X)
lines(t, beta[1,-1], col = "red")
lines(t, beta[1,-1] + 2*P[1,1,-1], lty = 2, col = "red")
lines(t, beta[1,-1] - 2*P[1,1,-1], lty = 2, col = "red")
lines(t, beta[1,-1] + 2 * (kappa[1,-1] + R[1,1,-1]), lty = 2, col = "blue")
lines(t, beta[1,-1] - 2 * (kappa[1,-1] + R[1,1,-1]), lty = 2, col = "blue")

