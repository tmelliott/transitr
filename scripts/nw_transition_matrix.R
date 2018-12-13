## create a transition matrix using adjacent segments
library(tidyverse)
library(ggraph)
library(tidygraph)
library(Matrix)

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

plotspeedmatrix <- function(x) {
    colnames(x) <- paste0("seg", 1:ncol(x))
    x %>% as.tibble %>% mutate(t = 1:n()) %>%
        gather(key = "segment", value = "speed", -t) %>%
        ggplot(aes(t, speed, group = segment, color = segment)) +
            geom_path()
}

### simulation 1: uncorrelated speeds
## cols: segments; rows: times
set.seed(12345)
Y <- matrix(NA, nrow = 100, ncol = 6)
Y[1, ] <- c(10, 8, 20, 7, 12, 14)
for (i in 2:nrow(Y)) {
    Y[i, ] <- rnorm(ncol(Y), Y[i-1,], 0.5)
}
plotspeedmatrix(Y)
apply(Y, 1, function(x) {
    dev.hold()
    print(drawnw(nwgraph, x))
    dev.flush(dev.flush())
    # Sys.sleep(1)
})

## we want to estimate p(F | mu, sigma, X, Y), assuming mu(t) = mu for all t
X <- Y * NA

jags.data <- list(M = ncol(Y), T = nrow(Y), Y = Y)

jags.fit <- jags.model('scripts/nw_model.jags', jags.data, n.chains = 4)
fsamps <- coda.samples(jags.fit, c('sigma'), n.iter = 10000, thin = 10)
summary(fsamps)

plot(fsamps)

F <- generate_matrix(segments)
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
