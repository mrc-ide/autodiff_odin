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

sim <- rstan::constrain_pars(stan_fit,
                      rstan::unconstrain_pars(stan_fit,list(beta=0.25, gamma = 0.1, I0=1)))

sir <- odin.dust::odin_dust("models/sir_4_AD.R")
beta <- 0.25
gamma <- 0.1
I0 <- 1
pars <- list(beta = beta, gamma = gamma, I0 = I0)
sir_model <- sir$new(pars, 0, 1)
y <- sir_model$simulate(seq(0,400, by=4))

plot(y[6,1,-1])
lines(sim$Y_model, col="red")

gradient_sir <- function(x, sir_gen, adj_gen, data_input){
  beta <- x[1]
  gamma <- x[2]
  I0 <- x[3]
  pars <- list(beta = beta, gamma = gamma, I0 = I0)
  sir_model <- sir_gen$new(pars, 0, 1)
  y <- sir_model$simulate(seq(0,400))
  pars <- list(beta = beta,
               gamma = gamma,
               I0 = I0,
               main_states = y[,1,],
               data_input = data_input,
               total_steps = 400)
  adj_model <- adj_gen$new(pars, 0, 1)
  adj_y <- adj_model$simulate(400)
  adj_y[7:9,1,1]
}

adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")

data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)
data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

gradient_sir(c(beta,gamma,I0), sir, adj_sir, data_input)
