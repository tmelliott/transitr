library(tidyverse)

## look at what the simulations look like
sim <- "sim000"
config <- jsonlite::read_json(file.path("simulations", sim, "config.json"))

transition <- function(state, delta, Sd = c(0, 1e5), config, id) {
    Dmax <- max(Sd)
    if (state[1] >= Dmax) {
        state[1] <- Dmax
        return(state)
    }
    
    M <- length(Sd)
    m <- which(Sd > state[1])[1] - 1
    next_stop_d <- Sd[m+1]

    State <- matrix(NA, nrow = delta+1, ncol = 4)
    while (state[1] < Dmax && delta > 0) {
        State[delta+1,] <- state

        n <- 0
        accel <- -100
        while ((state[2] + accel < 0 || state[2] + accel > 30)) {
            accel <- rnorm(1) * (1 + n/100) * config$system_noise

            n <- n + 1
        }
        
        # abar <- config$accel_rate * state[4]
        # alpha <- min(0, dnorm(accel, abar, 0.5, log = TRUE) -
        #     dnorm(state[3], abar, 0.5, log = TRUE))
        # if (runif(1) > exp(alpha)) {
        #     ## reject acceleration
        #     accel <- state[3]
        # }

        v <- max(0, min(30, state[2] + state[3]))
        vstar <- max(0, min(30, state[2] + accel))
        alpha <- min(0, dnorm(vstar, config$speed_mean, config$speed_sd, log = TRUE) - 
            dnorm(v, config$speed_mean, config$speed_sd, log = TRUE))
        if (runif(1) < exp(alpha)) {
            ## reject
            state[3] <- accel
            state[2] <- vstar
        } else {
            state[2] <- v
        }

        if (state[1] + state[2] >= next_stop_d) {
            if (runif(1) < 0.5) {
                dwell <- 6 + rexp(1, 1/10)
                delta <- max(0, delta - dwell)
                state[1] <- next_stop_d
                state[3] <- 0
            }
            m <- m+1
            if (m == M) break
            next_stop_d <- Sd[m+1]
            next
        }
        state[1] <- state[1] + state[2]
        delta <- delta - 1
    }
    State[1,] <- state
    colnames(State) <- c("distance", "speed", "acceleration", "state")
    State <- as.tibble(State[nrow(State):1,])
    if (!missing(id)) 
        State <- State %>% mutate(id = id)
    class(State) <- c("state.history", class(State))
    attr(state, "history") <- State
    class(state) <- "state"
    invisible(state)
}
plot.state <- function(x, ...) {
    plot(attr(x, "history"))
}
plot.state.history <- function(x, config, ...) {
    x <- x %>% group_by(id) %>% 
        do((.) %>% mutate(t = 1:n()-1)) %>%
        ungroup () %>% filter(!is.na(distance))
    p <- ggplot(x, aes(t/ifelse(max(x$t) < 60, 1, 60), group = id)) +
        xlab(ifelse(max(x$t) < 60, "Time (sec)", "Time (min)")) +
        theme(legend.position = "none")
    p1 <- p + geom_step(aes(y = distance))
    p2 <- p + geom_step(aes(y = speed))
    if (!missing(config)) {
        p2 <- p2 + 
            geom_hline(yintercept = c(config$speed_mean - 2 * config$speed_sd,
                                      config$speed_mean + 2 * config$speed_sd),
                        col = "orangered", lty = 3) +
            geom_hline(yintercept = c(config$speed_mean), col = "orangered")
    }
    p3 <- p + geom_step(aes(y = acceleration))
    p4 <- p + geom_step(aes(y = state))

    gridExtra::grid.arrange(p1, p2, ncol = 1)#, p3, p4, ncol = 1)
}

Stops <- c(0, 3000, 6000, 10000)
config$system_noise <- 0.5
config$accel_rate <- 1.5
config$speed_mean <- 13
config$speed_sd <- 2
x <- pbapply::pblapply(1:100, function(i) 
    attr(transition(c(0, 10, 0, 0), 600, Sd = Stops, config = config, id = i), "history")) %>%
    do.call(bind_rows, .)
plot(x, config)
