m <- rstan::stan_model("models/SIR.stan")

incidence <- read.csv("data/incidence.csv")

stan_fit <- rstan::sampling(m, data = list(
  T=100,
  Y=incidence$cases,
  freq=4),
  chains = 4, iter = 2000, thin = 1, cores = 4, control = list(adapt_delta = 0.95, max_treedepth = 15),
  init = list(list(beta=0.25, gamma = 0.1, I0=1),
              list(beta=1, gamma = 1, I0=10),
              list(beta=1, gamma = 1, I0=10),
              list(beta=1, gamma = 1, I0=10)))
multistrain_fit <- rstan::extract(stan_fit)

rstan::grad_log_prob(stan_fit, c(0.25, 0.1, 1))
