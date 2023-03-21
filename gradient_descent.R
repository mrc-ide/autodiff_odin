incidence <- read.csv("data/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

sir <- odin.dust::odin_dust("models/sir_4_AD.R")
adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")

compare <- function(state, observed, pars = NULL) {
  modelled <- state["incidence", , drop = TRUE]
  lambda <- modelled #+ rexp(length(modelled), 1e6)
  dpois(observed$cases, lambda, log = TRUE)
}

index <- function(info) {
  list(run = c(incidence = info$index$cases_inc),
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))
}

data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

gradient_sir <- function(beta = 0.25, gamma = 0.1, sir_gen, adj_gen, data_input){
  pars <- list(beta = beta, gamma = gamma, I0 = 1)
  sir_model <- sir_gen$new(pars, 0, 1)
  y <- sir_model$simulate(seq(0,400))
  y <- y[,1,]
  pars <- list(beta = beta,
               gamma = gamma,
               I0 = 1,
               main_states = y,
               data_input = data_input,
               total_steps = 400)
  adj_model <- adj_gen$new(pars, 0, 1)
  adj_y <- adj_model$simulate(seq(0,400))
  adj_y[7:8,1,401]
}

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
                                       compare = compare, index = index)

learning_rate <- 0.00002
beta_c <- 0.25
gamma_c <- 0.1
for(i in 1:20){
  g <- gradient_sir(beta = beta_c, gamma = gamma_c, sir, adj_sir, data_input)
  beta_c <- beta_c - g[1]*learning_rate
  gamma_c <- gamma_c - g[2]*learning_rate
  print(filter$run(list(beta = beta_c, gamma = gamma_c)))
}
