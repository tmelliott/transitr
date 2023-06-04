library(chopthin)
library(microbenchmark)

N <- 1e5
wt <- rbeta(N, 0.01, 1)
wt <- wt / sum(wt) * N


microbenchmark(
    sample(N, size = N, replace = TRUE, prob = wt),
    chopthin(wt, N)
)
