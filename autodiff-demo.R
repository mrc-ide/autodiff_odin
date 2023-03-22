incidence <- read.csv("data/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

#sir <- odin.dust::odin_dust("models/sir.R")

parsed_model <- jsonlite::fromJSON(odin::odin_parse("models/sir_4_AD.R"))

construct_param_tree <- function(parsed_model){
  parameter_graph <- NULL
  for(n in seq_along(parsed_model$equations$name)){
        if(!is.null(parsed_model$equations$depends$variables[n][[1]])){
          for(p in eval(parsed_model$equations$depends$variables[n][[1]]))
            parameter_graph <- rbind(parameter_graph,
                                     c(p,parsed_model$equations$name[n]))
          }
  }
  igraph::graph_from_edgelist(parameter_graph)
}

param_graph <- construct_param_tree(parsed_model)

plot.igraph(param_graph)

sir <- odin.dust::odin_dust("models/sir_4_AD.R")

pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
mod <- sir$new(pars, 0, 1)
y <- mod$simulate(seq(0,400))
i <- mod$info()$index[["time"]]
j <- mod$info()$index[["cases_inc"]]
plot(y[i, 1, ], t(y[j,1, ]), type = "l", col = "#00000055", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases", lwd=3)
points(cases ~ day, incidence, col = "red", pch = 19)



index <- function(info) {
  list(run = c(incidence = info$index$cases_inc),
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))
}

compare <- function(state, observed, pars = NULL) {
  dpois(observed$cases, state, log = TRUE)
}

LL <- 0

# compare <- function(state, observed, pars = NULL) {
#   modelled <- state["incidence", , drop = TRUE]
#   lambda <- modelled + rexp(length(modelled), 1e6)
#   dpois(observed$cases, lambda, log = TRUE)
# }

# filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 100,
#                                        compare = compare, index = index)
#
# priors <- list(
#   mcstate::pmcmc_parameter("beta", 0.2, min = 0),
#   mcstate::pmcmc_parameter("gamma", 0.1, min = 0, prior = function(p)
#     dgamma(p, shape = 1, scale = 0.2, log = TRUE)))
#
# vcv <- matrix(c(0.00057, 0.00052, 0.00052, 0.00057), 2, 2)
#
# transform <- function(theta) {
#   as.list(theta)
# }
#
# mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform)
#
# vcv <- matrix(c(0.00057, 0.00052, 0.00052, 0.00057), 2, 2)
# mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform)
# control <- mcstate::pmcmc_control(
#     n_steps = 500,
#     n_chains = 4,
#     n_threads_total = 12,
#     n_workers = 4,
#     save_state = TRUE,
#     save_trajectories = TRUE,
#     progress = TRUE)
# samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

sir <- odin.dust::odin_dust("models/sir_4_AD.R")
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
mod <- sir$new(pars, 0, 1)
y <- mod$simulate(seq(0,400))
y <- y[,1,]
adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")
pars <- list(beta = 0.25,
             gamma = 0.1,
             I0 = 1,
             main_states = y,
             data_input = data_input,
             total_steps = 401)
adj_mod <- adj_sir$new(pars, 0, 1)
adj_y <- adj_mod$simulate(seq(0,400))

data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

pars <- list(beta = 0.25, gamma = 0.1)
mod <- sir$new(pars, 0, 20)
y <- mod$simulate(c(0, data$time_end))
i <- mod$info()$index[["time"]]
j <- mod$info()$index[["cases_inc"]]
matplot(y[i, 1, ], t(y[j, , ]), type = "l", col = "#00000055", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases")
points(cases ~ day, incidence, col = "red", pch = 19)

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

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 100,
                                       compare = compare, index = index)
pars <- list(beta = 0.25, gamma = 0.1)
ff <- filter$run(pars)

pars <- list(beta = 0.25 + 1e-6, gamma = 0.1)
ff_beta <- filter$run(pars)

pars <- list(beta = 0.25, gamma = 0.1 + 1e-6)
ff_gamma <- filter$run(pars)

index(mod$info())
