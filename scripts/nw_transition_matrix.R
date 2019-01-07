## create a transition matrix using adjacent segments
library(tidyverse)
library(ggraph)
library(tidygraph)
library(Matrix)
library(rjags)
library(RSQLite)
library(dbplyr)


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
generate_matrix <- function(x, deg = 1) {
    m <- Matrix(diag(nrow(x)), doDiag = FALSE)
    if (deg == 0) return(m)
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
    diag(m) <- 0
    m <- sweep(m, 1, 2 * rowSums(m), "/")
    diag(m) <- 0.5
    m
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
Y <- X * NA
for (i in 1:nrow(X)) Y[i, ] <- rnorm(ncol(X), X[i, ], 1)
plotspeedmatrix(Y)

# 2. use data ("speeds") to calculate coefficients
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
fsamps <- coda.samples(jags.fit, 
    c('alpha_raw', 'alpha', 'pr_alpha_zero', 'sigma', 'q'), 
    n.iter = 1000, thin = 1)
summary(fsamps)

devAskNewPage(TRUE)
plot(fsamps)
devAskNewPage(FALSE)

## extract a simulation and plot
sm <- coda.samples(jags.fit, c('X'), n.iter = 1000, thin = 10) %>% 
    as.matrix
for (i in 1:100) {
    dev.hold()
    print(plotspeedmatrix(Y, sm[sample(1:nrow(sm), 1), ]))
    dev.flush()
}

Xsamps <- coda.samples(jags.fit, c('alpha'), n.iter = 10000, thin = 10)
plot(Xsamps)


### Use some real data!

get_segment_data <- function(routes) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    if (missing(routes)) {
        segments <- con %>% tbl("road_segments")
    } else {
        rids <- con %>% tbl("routes") %>%
            filter(route_short_name %in% routes &
                (route_long_name %like% "%To City%" |
                 route_long_name %like% "%To Britomart%" |
                 route_long_name %like% "%To Mayoral%" |
                 route_long_name %like% "%To Auckland Universities%")) %>%
            collect %>% pluck("route_id") %>% unique
        sids <- con %>% tbl("trips") %>%
            filter(route_id %in% rids) %>%
            collect %>% pluck("shape_id") %>% unique
        segids <- con %>% tbl("shape_segments") %>% 
            filter(shape_id %in% sids) %>%
            collect %>% pluck("road_segment_id") %>% unique
        segments <- con %>% tbl("road_segments") %>%
            filter(road_segment_id %in% segids)
    }
    intersections <- con %>% tbl("intersections")
    segments <- segments %>% 
        inner_join(intersections, by = c("int_from" = "intersection_id"), suffix = c("", "_start")) %>%
        inner_join(intersections, by = c("int_to" = "intersection_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, intersection_lat, intersection_lon, intersection_lat_end, intersection_lon_end) %>%
        collect
    segments
}

read_segment_data <- function(sim) {
    file.path("simulations", sim, "history") %>%
        list.files(pattern = "segment_", full.names = TRUE) %>%
        lapply(
            read_csv, 
            col_types = "ciid", 
            col_names = 
                c("segment_id", "timestamp", "travel_time", "error")
        ) %>%
        bind_rows() %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))
}

## just segments of these routes
segdata <- get_segment_data(c("27W", "27H", "22N", "22R", "22A"))
segdata <- get_segment_data(c("NX1", "NX2", "82"))
sids <- unique(segdata$road_segment_id)

con <- dbConnect(SQLite(), "fulldata.db")
seglens <- con %>% tbl("road_segments") %>% 
    filter(road_segment_id %in% sids) %>%
    select(road_segment_id, length, int_from, int_to) %>% collect %>%
    mutate(road_segment_id = as.character(road_segment_id),
           from = int_from, to = int_to)
dbDisconnect(con)


rawdata <- read_segment_data("sim000") %>% 
    filter(segment_id %in% sids) %>%
    left_join(seglens, by = c("segment_id" = "road_segment_id"))
segg <- rawdata %>%
    mutate(speed = length / travel_time) %>%
    group_by(segment_id) %>% 
        summarize(
            speed = mean(speed) * 60 * 60 / 1000,
            from = first(from), to = first(to)
        ) %>%
    as_tbl_graph()

ggraph(segg, layout = 'kk') + 
    geom_edge_fan(
        aes(color = speed),
        arrow = arrow(length = unit(1, 'mm')), 
        start_cap = circle(1, 'mm'),
        end_cap = circle(1, 'mm')) +
    scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limits = c(0, 100)) +
    geom_node_point(shape = 19, size = 0.5)

### by the hour (for now)
discretise <- function(x, seconds = 3600) {
    date <- as.integer(as.POSIXct(format(x$timestamp[1], "%Y-%m-%d 00:00:00")))
    bins <- seq(5*60*60, 24*60*60, by = seconds)
    ts <- as.integer(x$timestamp) - date
    x$time <- cut(ts, breaks = bins, 
        labels = as.POSIXct(date + bins[-length(bins)], origin = "1970-01-01") %>%
            format("%T"))
    x
}

seg.hour <- rawdata %>% discretise(60*15) %>%
    mutate(speed = length / travel_time) %>%
    group_by(segment_id, time) %>% 
        summarize(
            speed = mean(speed) * 60 * 60 / 1000,
            from = first(from), to = first(to)
        )
segg.hour <- seg.hour %>% as_tbl_graph() %>% activate(edges)

ggraph(segg.hour, layout = 'kk') + 
    geom_edge_fan(aes(color = speed)) +
    scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limits = c(0, 100)) +
    facet_wrap(~time)

## Generate observation matrix
segY <- seg.hour %>% spread(key = time, value = speed) %>% ungroup() %>%
    mutate(segment_id = as.character(segment_id)) %>% arrange(segment_id)
Y <- segY %>% select(-from, -to, -segment_id) %>% as.matrix

segm <- segY %>% mutate(id = segment_id) %>% select(id, from, to)

## independent version
F <- generate_matrix(segm, 0) %>% diag()
jags.data <- list(
    M = nrow(Y),
    T = ncol(Y),
    Y = Y,
    alpha = F
)




## non-independence version
F <- generate_matrix(segm)

tF <- t(F) ## flip to get row-major form, for easier manipulation 
jags.data <- list(
    M = nrow(Y), 
    T = ncol(Y), 
    Y = Y,
    alpha = attr(tF, "x"),
    Fp = attr(tF, "p"),     # the column indexes
    Fi = attr(tF, "i") + 1 # the row indexes
    # NNZ = length(attr(tF, "x"))
)

jags.fit <- jags.model('scripts/nw_model.jags', jags.data, 
    n.chains = 1, n.adapt = 10000)
fsamps <- coda.samples(jags.fit, 
    c('sigma', 'q'), 
    # c('alpha_raw', 'pr_alpha_zero', 'sigma', 'q'), 
    n.iter = 10000, thin = 10)
summary(fsamps)

plot(fsamps)

betas <- coda.samples(jags.fit, "beta", n.iter = 10000, thin = 10)
Fhat <- tF
attr(Fhat, "x") <- as.matrix(betas) %>% apply(1, median)
Fhat <- t(Fhat)
image(Fhat, cuts = 20, col.regions = viridis::viridis(21), 
    colorkey = TRUE, lwd = 0)



mb <- as.matrix(betas)
tfhat <- t(F)
for (i in sample(1:nrow(mb), 100)) {
    attr(tfhat, "x") <- mb[i, ]
    Fhat <- t(tfhat)
    dev.hold()
    print(image(Fhat, cuts = 20, col.regions = viridis::viridis(21), 
        colorkey = TRUE, lwd = 0))
    dev.flush()
}

devAskNewPage(TRUE)
plot(betas)
devAskNewPage(FALSE)

Xsamps <- coda.samples(jags.fit, "X", n.iter = 1000)

drawX <- function(x, d, t = "12:00:00") {
    ## put x into a matrix
    dx <- expand.grid(
        segment_id = unique(d$segment_id), 
        time = unique(d$time)
    ) %>% 
        as.tibble %>%
        mutate(speed = x, segment_id = as.character(segment_id))
    ds <- d %>% summarize(from = first(from), to = first(to))
    g <- dx %>% left_join(ds) %>%
        as_tbl_graph %>% activate(edges)
    ggraph(g %>% filter(time == t), layout = 'kk') + 
        geom_edge_fan(
            aes(color = speed),
            arrow = arrow(length = unit(1, 'mm')), 
            edge_width = 1,
            start_cap = circle(1, 'mm'),
            end_cap = circle(1, 'mm')) +
        scale_edge_colour_gradientn(
            colours = viridis::viridis(256, option = "A"), 
            limits = c(0, 100)
        ) +
        geom_node_point(shape = 19, size = 0.5)
}

Xm <- as.matrix(Xsamps)
dev.flush(dev.flush())
si <- sample(1:nrow(Xm), 1)
for (t in sort(unique(seg.hour$time))) {
    dev.hold()
    print(drawX(Xm[si, ], seg.hour, t) + ggtitle(t))
    dev.flush()
}

for (i in sample(1:nrow(Xm), 100)) {
    Xmi <- expand.grid(
            segment_id = sort(unique(seg.hour$segment_id)),
            time = sort(unique(seg.hour$time))
        ) %>%
        as.tibble %>%
        mutate(segment_id = as.character(segment_id)) %>%
        left_join(seg.hour %>% summarize(from = first(from), to = first(to))) %>%
        mutate(speed = Xm[i, ],
            time = as.POSIXct(paste(Sys.Date(), time), origin = "1970-01-01"))
    dev.hold()
    print(ggplot(Xmi, aes(time, speed, group = segment_id)) +
        geom_path(aes(colour = segment_id)) +
        ylim(0, 100))
    dev.flush(dev.flush())
}
