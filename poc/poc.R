gen <- dust::dust("adjoint.cpp")

## The first thing I need to do is to get the likelihood calculation
## working at all from dust. We should be able to do this with the
## filter.
incidence <- read.csv("marc/incidence.csv")
d <- dust::dust_data(
  data.frame(incidence = incidence$cases, time = incidence$day * 4))
m1 <- gen$new(list(beta = 0.25, gamma = 0.1, I0 = 1), 0, 1,
              deterministic = TRUE)
m1$set_data(d)
m1$filter()$log_likelihood

log_likelihood <- function(..., m) {
  m$update_state(pars = list(...), time = 0)
  m$filter()$log_likelihood
}

beta <- 0.25
gamma <- 0.1
I0 <- 1

log_likelihood(beta = beta, gamma = gamma, I0 = I0, m = m1) # -482.4857

h <- 1e-6
(log_likelihood(beta = beta + h, gamma = gamma, I0 = I0, m = m1) -
 log_likelihood(beta = beta - h, gamma = gamma, I0 = I0, m = m1)) / (2 * h)
(log_likelihood(beta = beta, gamma = gamma + h, I0 = I0, m = m1) -
 log_likelihood(beta = beta, gamma = gamma - h, I0 = I0, m = m1)) / (2 * h)

(log_likelihood(beta = beta + h, gamma = gamma, I0 = I0, m = m1) -
 log_likelihood(beta = beta, gamma = gamma, I0 = I0, m = m1)) / h
(log_likelihood(beta = beta, gamma = gamma + h, I0 = I0, m = m1) -
 log_likelihood(beta = beta, gamma = gamma, I0 = I0, m = m1)) / h

## This runs ok, but produces the wrong values. This might be because:
## * I could have implemented the update function wrong
## * I could have off-by-ones (or misses) in passing the data in
##
## now that we have the above basically working, compute a 4 step, one
## data point case and work back from that.
res <- gen$parent_env$newthing(list(beta = beta, gamma = gamma, I0 = I0), d)
res[[2]][6:8]

m2 <- gen$new(list(beta = beta, gamma = gamma, I0 = I0), 0, 1,
              deterministic = TRUE)
m2$set_data(data)
bench::mark(
  gen$parent_env$newthing(m2$.__enclos_env__$private$ptr_),
  log_likelihood(m = m),
  check = FALSE)




## ------


pkgload::load_all(compile = FALSE, recompile = FALSE)
gen <- dust(dust_file("examples/ode/logistic.cpp"))
pars <- list(r = c(0.1, 0.2), K = c(100, 100))
control <- dust_ode_control(step_size_min = 0.1)
n_particles <- 5
mod <- gen$new(pars, 0, n_particles, ode_control = control)
expect_error(
  mod$run(5),
  "step too small")
expect_error(
  mod$run(5),
  "Errors pending; reset required")
mod$update_state(pars = pars, reset_step_size = TRUE)


gen <- dust(dust_file("examples/sirs.cpp"))

## pkgload::load_all(compile = FALSE, recompile = FALSE)
gen <- dust("inst/examples/ode/malaria.cpp")
mod <- gen$new(list(), 0, 50, seed = 1L)
cases <- head(read.csv("inst/extdata/malaria_cases.csv"), 40)
mod$set_stochastic_schedule(cases$t[-1])
y <- mod$simulate(cases$t)
#$ saveRDS(y, "y.rds")
cmp <- readRDS("y.rds")
all.equal(y, cmp)


pkgload::load_all(compile = FALSE, recompile = FALSE)
gen <- dust(dust_file("examples/sirs.cpp"), skip_cache = TRUE)


gen <- odin::odin({
  initial(x) <- 1
  update(x) <- sum(y)
  y[f(1:n + 1)] <- 1
  dim(y) <- m
  n <- 5
  m <- n + 1
})

