gen <- dust::dust("adjoint.cpp")

## The first thing I need to do is to get the likelihood calculation
## working at all from dust. We should be able to do this with the
## filter.
read_data <- function(n = 100, freq = 4) {
  incidence <- read.csv("marc/incidence.csv")
  incidence <- incidence[incidence$day <= n, ]
  dust::dust_data(
    data.frame(incidence = incidence$cases, time = incidence$day * freq))
}

d <- read_data()

pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
## m1 <- gen$new(pars, 0, 1, deterministic = TRUE)
## m1$set_data(d)
## m1$filter()$log_likelihood

res <- gen$parent_env$newthing(pars, d)
res

expected <- c(4511.91973090839, -2892.37818928496, 234.49219517371)
expected[[3]] <- 0
testthat::expect_equal(res[[2]], expected, tolerance = 1e-12)
