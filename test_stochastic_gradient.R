incidence <- read.csv("data/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

sir_stoch <- odin.dust::odin_dust("models/sir_stochastic.R")

sir_stoch_model <- sir_stoch$new(pars, 0, 100)
beta <- 0.25
gamma <- 0.1
I0 <- 1
pars <- list(beta = beta, gamma = gamma, I0 = I0)
