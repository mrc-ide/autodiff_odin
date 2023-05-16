## This file pulls out the core bits of the qmd file into an easily
## runnable script so that we can then establish equivalence with the
## proof of concept.
sir <- odin.dust::odin_dust("marc/sir.R")
sir_adjoint <- odin.dust::odin_dust("marc/sir_adjoint.R")

index <- function(info) {
  list(run = c(incidence = info$index$cases_inc),
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))
}

compare <- function(state, observed, pars = NULL) {
  modelled <- state["incidence", , drop = TRUE]
  lambda <- modelled
  dpois(observed$cases, lambda, log = TRUE)
}

beta <- 0.25
gamma <- 0.1
I0 <- 1
pars <- list(beta = beta, gamma = gamma, I0 = I0)
sir_model <- sir$new(pars, 0, 1)
y <- sir_model$simulate(seq(0, 400))

incidence <- read.csv("marc/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

data_input <- NULL
for (i in 1:100) {
  data_input <- c(data_input, rep(0,3), data$cases[i])
}

y <- y[, 1, ]
pars <- list(beta = beta,
             gamma = gamma,
             I0 = I0,
             main_states = y,
             data_input = data_input,
             total_steps = 400)
adj_model <- sir_adjoint$new(pars, 0, 1)
adj_y <- adj_model$simulate(400)
adj_y[7:9, 1, 1]

gradient_sir <- function(x, sir_gen, adj_gen, data_input) {
  beta <- x[1]
  gamma <- x[2]
  I0 <- x[3]
  pars <- list(beta = beta, gamma = gamma, I0 = I0)
  sir_model <- sir_gen$new(pars, 0, 1)
  y <- sir_model$simulate(seq(0, 400))
  pars <- list(beta = beta,
               gamma = gamma,
               I0 = I0,
               main_states = y[,1,],
               data_input = data_input,
               total_steps = 400)
  adj_model <- adj_gen$new(pars, 0, 1)
  adj_y <- adj_model$simulate(400)
  adj_y[7:9, 1, 1]
}

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
                                       compare = compare, index = index)

num_grad_sir <- function(x, filter, h = 1e-6) {
  beta <- x[1]
  gamma <- x[2]
  I0 <- x[3]
  ll_c <- filter$run(list(beta = beta, gamma = gamma, I0 = I0))
  ll_beta_h <- filter$run(list(beta = beta+h, gamma = gamma, I0 = I0))
  ll_gamma_h <- filter$run(list(beta = beta, gamma = gamma+h, I0 = I0))
  ll_I0_h <- filter$run(list(beta = beta, gamma = gamma, I0 = I0+h))
  list(log_likelihood = ll_c,
       gradient = (c(ll_beta_h, ll_gamma_h, ll_I0_h) - ll_c) / h)
}

num_grad_sir(c(beta, gamma, I0), filter)
gradient_sir(c(beta, gamma, I0), sir, sir_adjoint, data_input)

x_c_0 <- c(.35, 0.2, 8)
bench::mark(
  num_grad_sir(x_c_0, filter),
  gradient_sir(x_c_0, sir, sir_adjoint, data_input),
  check = FALSE)

## with beta = 0.25, gamma = 0.1, I0 = 1 -> ll = -482.4857
m <- sir$new(list(beta = beta, gamma = gamma, I0 = I0), 0, 1,
             deterministic = TRUE)
d <- dust::dust_data(
  data.frame(cases_inc = incidence$cases, time = incidence$day * 4))
m$set_data(d)
m$filter()$log_likelihood
