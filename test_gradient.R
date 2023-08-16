gen <- odin.dust::odin_dust("models/sir_adjoint.R")

incidence <- read.csv("data/incidence.csv")
incidence <- data.frame(
  time = incidence$day * 4,
  cases_observed = incidence$cases)
d <- dust::dust_data(incidence)

pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
mod <- gen$new(pars, 0, 1, deterministic = TRUE)
mod$set_data(d)

## This is the current temporary arrangement with dust and may change:
info <- mod$info()


plot(incidence)

mod$update_state(list(beta = 0.7, gamma = 0.1, I0 = 1), time = 0)
mod$set_data(d)
mod$run_adjoint()

x <- c(0.25,0.1,1)
h <- 1e-7
names(x) <- c("beta","gamma","I0")
LL_chain <- NULL

for(i in 1:1000){
  res <- mod$run_adjoint()
  LL_chain <- c(LL_chain, res$log_likelihood)
  x <- x + h * res$gradient
  mod$update_state(list(beta = x["beta"], gamma = x["gamma"], I0 = x["I0"]), time = 0)
}

x <- c(0.25, 0.1, 1)
h <- 1e-8
names(x) <- c("beta", "gamma", "I0")
LL_chain2 <- NULL

# Initialize NAG-specific variables
momentum <- 0.8
velocity <- rep(0, length(x))

for (i in 1:1000) {
  # Calculate gradient using the current lookahead position

  lookahead_x <- x + momentum * velocity
  mod$update_state(list(beta = lookahead_x["beta"], gamma = lookahead_x["gamma"], I0 = lookahead_x["I0"]), time = 0)
  res <- mod$run_adjoint()
  gradient <- res$gradient

  # Update Nesterov velocity
  velocity <- momentum * velocity + h * gradient

  # Update the parameters using the Nesterov velocity
  x <- x + velocity

  # Update model state
  mod$update_state(list(beta = x["beta"], gamma = x["gamma"], I0 = x["I0"]), time = 0)
  res <- mod$run_adjoint()
  LL_chain2 <- c(LL_chain2, res$log_likelihood)
}
