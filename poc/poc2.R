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

d <- read_data(10)
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
res <- gen$parent_env$newthing(pars, d)
res

expected10 <- list(-44.0256051296862,
                   c(244.877646917118, -140.566517375877,  25.2152128116894))
testthat::expect_equal(res, expected10, tolerance = 1e-12)



expected <- c(4511.91973090839, -2892.37818928496, 234.49219517371)
testthat::expect_equal(res[[2]], expected, tolerance = 1e-12)
