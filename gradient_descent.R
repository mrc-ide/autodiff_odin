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

gradient_sir <- function(x, sir_gen, adj_gen, data_input){
  beta <- x[1]
  gamma <- x[2]
  I0 <- x[3]
  pars <- list(beta = beta, gamma = gamma, I0 = I0)
  sir_model <- sir_gen$new(pars, 0, 1)
  y <- sir_model$simulate(seq(0,400))
  y <- y[,1,]
  pars <- list(beta = beta,
               gamma = gamma,
               I0 = I0,
               main_states = y,
               data_input = data_input,
               total_steps = 400)
  adj_model <- adj_gen$new(pars, 0, 1)
  adj_y <- adj_model$simulate(seq(0,400))
  adj_y[7:9,1,401]
}

num_grad_sir <- function(x, filter, h = 1e-6){
  beta <- x[1]
  gamma <- x[2]
  I0 <- x[3]
  LL_c <- filter$run(list(beta = beta, gamma = gamma, I0 = I0))
  LL_beta_h <- filter$run(list(beta = beta+h, gamma = gamma, I0 = I0))
  LL_gamma_h <- filter$run(list(beta = beta, gamma = gamma+h, I0 = I0))
  LL_I0_h <- filter$run(list(beta = beta, gamma = gamma, I0 = I0+h))
  (c(LL_beta_h,LL_gamma_h, LL_I0_h)-LL_c)/h
}

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
                                       compare = compare, index = index)
chain_values <- NULL
learning_rate <- 0.000008

x_c <- c(0.25,0.1,10)
for(i in 1:100){
  g <- gradient_sir(x_c, sir, adj_sir, data_input)
  #g_n <- num_grad_sir(x_c, filter, h = 1e-6)
  x_c <- x_c + g*learning_rate
  LL <- filter$run(list(beta = x_c[1], gamma = x_c[2], I0 = x_c[3]))
  chain_values <- rbind(chain_values, c(x_c, LL))
}

plot(chain_values[,4], type="l")

