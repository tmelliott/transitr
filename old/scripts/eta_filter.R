library(tidyverse)
## a (constrained) KF to ensure increasing ETAs along stops

F <- function(B, delta = 0) {
    x <- diag(nrow(B))
    if (delta == 0) return(x)
    for (i in 1:nrow(B)) {
        if (x[i, i] < 1e-6) continue()
        x[i, i] <- min(1, max(0, 1 - delta / B[i, ]));
    }
    x
}
Q <- function(B, delta) diag(nrow(B)) * delta / 2
H <- function(B) {
    x <- diag(nrow(B))
    for (i in 2:nrow(B)) x[i,i-1] = -1
    x
}

B0 <- cbind(1:5 * 30)
P0 <- diag(1:5 * 10)

data <- 
    list(
        list(
            delta = 20,
            Z = cbind(c(20, 55, 88, 115, 140)),
            R = diag(c(3, 6, 12, 20, 30))
        ),
        list(
            delta = 20,
            Z = cbind(c(0, 20, 55, 90, 115)),
            R = diag(c(0, 2, 10, 11, 15))
        ),
        list(
            delta = 20,
            Z = cbind(c(0, 0, 20, 55, 70)),
            R = diag(c(0, 0, 2, 3, 4))
        ),
        list(
            delta = 10,
            Z = cbind(c(0, 0, 14, 50, 60)),
            R = diag(c(0, 0, 2, 4, 8))
        ),
        list(
            delta = 10,
            Z = cbind(c(0, 0, 10, 40, 60)),
            R = diag(c(0, 0, 4, 5, 10))
        ),
        list(
            delta = 10,
            Z = cbind(c(0, 0, 2, 34, 52)),
            R = diag(c(0, 0, 1, 4, 12))
        )
    )

B <- B0
P <- P0

hist <- list()
t0 <- 0
for (k in 1:length(data)) {
    delta <- data[[k]]$delta
    Z <- data[[k]]$Z
    R <- data[[k]]$R^2

    Fk <- F(B, delta)
    Qk <- F(B, delta)
    Hk <- F(B, delta)
    
    ## kf prediction
    Bhat <- Fk %*% B
    Phat <- Fk %*% P %*% t(Fk) + Qk
    plot(B, 1:nrow(B), pch = 19, col = "gray", xlim = c(0, max(B0)), ylim = c(0, nrow(B)))
    arrows(B - sqrt(diag(P)), 1:nrow(B), B + sqrt(diag(P)), code = 0, col = "gray", lwd = 2)
    points(Z, 1:nrow(B) - 0.1, pch = 19, cex = 0.5, col = "black")
    arrows(Z - sqrt(diag(R)), 1:nrow(B) - 0.1, Z + sqrt(diag(R)), code = 0, col = "black")
    points(Bhat, 1:nrow(B) - 0.2, pch = 19, col = "red", cex = 0.8)
    arrows(Bhat - sqrt(diag(Phat)), 1:nrow(B) - 0.2, Bhat + sqrt(diag(Phat)), code = 0, col = "red")

    ## kf update
    yk <- Z - Hk %*% Bhat
    Sk <- R + Hk %*% Phat %*% t(Hk)
    Ski <- which(Z > 0)
    Kk <- Phat[Ski,Ski] %*% t(Hk[Ski,Ski]) %*% solve(Sk[Ski,Ski])
    B[-Ski] <- 0
    B[Ski] <- Bhat[Ski] + Kk %*% yk[Ski]
    P <- P * 0
    P[Ski,Ski] <- (1 - Kk %*% Hk[Ski,Ski]) %*% Phat[Ski,Ski] %*% t(1 - Kk %*% Hk[Ski,Ski]) + Kk %*% R[Ski,Ski] %*% t(Kk)
    points(B, 1:nrow(B) - 0.3, pch = 19, col = "blue", cex = 0.8)
    arrows(B - sqrt(diag(P)), 1:nrow(B) - 0.3, B + sqrt(diag(P)), code = 0, col = "blue")
    
    t0 <- t0 + delta
    hist[[k]] <- tibble(
        stop = 1:nrow(Z),
        time = t0,
        obs = Z %>% as.numeric,
        err = diag(R),
        eta_pred = Hk %*% Bhat %>% as.numeric,
        var_pred = diag(Hk %*% Phat %*% t(Hk)),
        eta = Hk %*% B %>% as.numeric,
        var = diag(Hk %*% Phat %*% t(Hk))
    )
}

## into tibble format
Hist <- bind_rows(hist)

ggplot(Hist, aes(stop)) +
    facet_grid(time~.) +
    geom_path(aes(y = obs), colour = "black") +
    geom_pointrange(aes(
        y = obs, ymin = obs - sqrt(err), ymax = obs + sqrt(err)
    ), pch = 21, fill = "gray") +
    geom_path(aes(y = eta_pred), colour = "red") +
    geom_pointrange(aes(
        y = eta_pred, ymin = eta_pred - sqrt(var_pred), ymax = eta_pred + sqrt(var_pred)
    ), colour = "red") +
    geom_path(aes(y = eta), colour = "blue") +
    geom_pointrange(aes(
        y = eta, ymin = eta - sqrt(var), ymax = eta + sqrt(var)
    ), colour = "blue")
