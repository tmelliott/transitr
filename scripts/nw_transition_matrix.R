## create a transition matrix using adjacent segments
library(tidyverse)
library(ggraph)
library(tidygraph)
library(Matrix)
library(rjags)

# 1. create a toy network
segments <- 
    tibble(
        id = 1:6,
        from = c(1:3, 5, 3, 6),
        to = c(2:4, 2, 6, 7)
    )

# 1b. plot the example
nwgraph <- as_tbl_graph(segments)
drawnw <- function(nw, speeds = 10, .seed = 1234) {
    set.seed(.seed)
    ggraph(nwgraph, layout="nicely") +
        geom_edge_fan(
            aes(color = speeds),
            arrow = arrow(length = unit(2, 'mm')), 
            start_cap = circle(2, 'mm'),
            end_cap = circle(2, 'mm')
        ) +
        geom_node_point() +
        scale_edge_color_gradientn(colours = viridis::viridis(256), limits = c(0, 30)) +
        theme(legend.position = "none")
}
drawnw(nwgraph)

# 2. generate matrix
generate_matrix <- function(x) {
    m <- Matrix(diag(nrow(x)), doDiag = FALSE)
    ## convert id to factor
    x$id <- as.integer(as.factor(x$id))
    ids <- unique(c(x$from, x$to))
    x$from <- as.integer(factor(x$from, levels = ids))
    x$to <- as.integer(factor(x$to, levels = ids))
    for (i in seq_along(1:nrow(m))) {
        toi <- x$id[x$to == x$from[i]]
        fromi <- x$id[x$from == x$to[i]]
        if (length(toi) > 0) m[i, toi] <- 1
        if (length(fromi) > 0) m[i, fromi] <- 1
    }
    sweep(m, 1, rowSums(m), "/")
}

mat <- generate_matrix(segments)

# 2b. test it
library(mvtnorm)

plotspeedmatrix <- function(x, fit) {
    colnames(x) <- paste0("seg", 1:ncol(x))
    p <- x %>% as.tibble %>% mutate(t = 1:n()) %>%
        gather(key = "segment", value = "speed", -t) %>%
        ggplot(aes(t, speed, group = segment, color = segment)) +
        ylim(range(x))
    if (missing(fit)) return(p + geom_path())

    ## fit should be a vector of values
    Xhat <- fit[grep("^X", names(fit))]
    Xhat <- t(matrix(Xhat, ncol = nrow(x)))
    colnames(Xhat) <- colnames(x)
    Xdf <- Xhat %>% as.tibble %>% mutate(t = 1:n()) %>%
        gather(key = "segment", value = "speed", -t)
    p + geom_path(lty = 2) +
        geom_path(data = Xdf)
}

### simulation 1: uncorrelated speeds
## cols: segments; rows: times
set.seed(12345)
X <- matrix(NA, nrow = 100, ncol = 6)
Xdot <- X
X[1, ] <- c(10, 8, 20, 7, 12, 14)
Xdot[1, ]<- c(0.01, -0.05, -1, 0.05, 0, -0.02)/10
F <- generate_matrix(segments)
tF <- t(F)
attr(tF, "x") <- c(
    1, 0, 
    0, 0.95, 0, 0.03, 0.02,
    0.03, 0.97,
    0, 1,
    0.06, 0.9, 0.04,
    0.01, 0.99
)
F <- t(tF)
for (i in 2:nrow(X)) {
    Xdot[i, ] <- rnorm(ncol(X), drop(F %*% Xdot[i-1,]), 0.002)
    X[i, ] <- X[i-1, ] + Xdot[i, ]
}
for (i in 1:nrow(X)) Y[i, ] <- rnorm(ncol(X), X[i, ], 1)
plotspeedmatrix(Y)
# apply(Y, 1, function(x) {
#     dev.hold()
#     print(drawnw(nwgraph, x))
#     dev.flush(dev.flush())
#     # Sys.sleep(1)
# })

## we want to estimate p(F | mu, sigma, X, Y), assuming mu(t) = mu for all t
X <- Y * NA
F <- generate_matrix(segments)
tF <- t(F) ## flip to get row-major form, for easier manipulation 

jags.data <- list(
    M = ncol(Y), 
    T = nrow(Y), 
    Y = t(Y), ## use transpose 
    # Fx = attr(tF, "x"),
    Fp = attr(tF, "p"),     # the column indexes
    Fi = attr(tF, "i") + 1, # the row indexes
    NNZ = length(attr(tF, "x"))
)

jags.fit <- jags.model('scripts/nw_model.jags', jags.data, 
    n.chains = 2, n.adapt = 10000)
fsamps <- coda.samples(jags.fit, c('X'), n.iter = 10000, thin = 10)
summary(fsamps)

asamps <- coda.samples(jags.fit, c('alpha'), n.iter = 10000, thin = 10)
plot(asamps)

# o <- par(mfrow = c(4, 6))
# traceplot(fsamps)
# par(o)

## extract a simulation and plot
sm <- as.matrix(fsamps)
for (i in 1:100) {
    dev.hold()
    print(plotspeedmatrix(Y, sm[sample(1:nrow(sm), 1), ]))
    dev.flush()
}




Q <- Matrix(diag(ncol(X)) * 10)  ## system noise
R <- Matrix(diag(ncol(X)) * 5)   ## measurement error
I <- Matrix(diag(ncol(X)))       ## an identity matrix




xhat <- cbind(rep(10, ncol(X)))
P <- Matrix(diag(ncol(X)) * 100, doDiag = FALSE)

pb <- txtProgressBar(0, nrow(Y), style = 3)
for (i in 1:nrow(Y)) {
    ## predict
    xhat <- as.matrix(F %*% xhat)
    P <- F %*% P %*% t(F) + Q

    ## update
    y <- Y[i, ] - xhat
    y[is.na(y)] <- 0
    S <- R + P
    K <- P %*% solve(S)
    X[i, ] <- xhat <- xhat + as.matrix(K %*% y)
    P <- (I - K) %*% P %*% t(I - K) + K %*% R %*% t(K)
    setTxtProgressBar(pb, i)
}; close(pb)


drawnw(nwgraph, speeds[, 1])
for (i in 1:nrow(X)) {
    p <- drawnw(nwgraph, X[i, ])
    dev.hold()
    print(p)
    dev.flush()
}

plotspeedmatrix(Y)


# 3. use data ("speeds") to calculate coefficients
