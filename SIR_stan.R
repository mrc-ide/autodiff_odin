m <- rstan::stan_model("models/SIR.stan")

incidence <- read.csv("data/incidence.csv")

n_data <- 100
stan_fit <- rstan::sampling(m, data = list(
  T=n_data,
  Y=incidence$cases[seq(n_data)],
  freq=4),
  chains = 4, iter = 2000, thin = 1, cores = 4, control = list(adapt_delta = 0.95, max_treedepth = 15),
  init = list(list(beta=0.25, gamma = 0.1, I0=1),
              list(beta=1, gamma = 1, I0=10),
              list(beta=1, gamma = 1, I0=10),
              list(beta=1, gamma = 1, I0=10)))
SIR_fit <- rstan::extract(stan_fit)

rstan::grad_log_prob(stan_fit,
                     rstan::unconstrain_pars(stan_fit,list(beta=0.25, gamma = 0.1, I0=1)))
