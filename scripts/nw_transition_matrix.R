## create a transition matrix using adjacent segments
library(tidyverse)
library(ggraph)
library(tidygraph)

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
        theme(legend.position = "none")
}
drawnw(nwgraph)

# 2. generate matrix
generate_matrix <- function(x) {
    m <- diag(nrow(x))
    for (i in seq_along(1:nrow(m))) {
        toi <- x$id[x$to == x$from[i]]
        fromi <- x$id[x$from == x$to[i]]
        m[i, toi] <- 1
        m[i, fromi] <- 1
    }
    sweep(m, 1, rowSums(m), "/")
}
mat <- generate_matrix(segments)

# 2b. test it
speeds <- cbind(c(10, 8, 20, 7, 12, 14))

drawnw(nwgraph, speeds)
for (i in 1:100) {
    p <- drawnw(nwgraph, speeds <- mat %*% speeds)
    dev.hold()
    print(p)
    dev.flush()
}



# 3. use data ("speeds") to calculate coefficients