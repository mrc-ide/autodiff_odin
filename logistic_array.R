logistic <- odin2::odin({
  update(x[]) <- r[i] * x[i] * (1 - x[i])

  initial(x[]) <- x0[i]

  x0 <- parameter()
  r <- parameter()
  N <- parameter()

  dim(x) <- 2
  dim(x0) <- 2
  dim(r) <- 2

  y <- sum(x)
  pop <- data()
  pop ~ Poisson(y * N)
})

pars <- list(
  x0 = c(0.001,.1),
  r = c(1.1, 1.3),
  N = 1000)

sys <- dust2::dust_system_create(sir, pars, dt = 1)

dust2::dust_system_set_state_initial(sys)
dust2::dust_system_state(sys)
t <- seq(1, 100)
y <- dust2::dust_system_simulate(sys, t)

set.seed(42)
pop <- rpois(length(t), pars$N * colSums(y[,seq_along(t)]))
data <- data.frame(time = t, pop = pop)

# Filter and likelihood
filter <- dust2::dust_unfilter_create(logistic, 0, data)

dust2::dust_likelihood_run(filter, pars)

plot(data$pop, pch=19, col= "red")
lines(colSums(y) * pars$N)

