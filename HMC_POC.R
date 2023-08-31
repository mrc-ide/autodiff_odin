# remotes::install_github("mrc-ide/odin@mrc-4358", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/odin.dust@mrc-4359", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/dust@mrc-4307", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/mcstate", upgrade = FALSE, force = TRUE)

# Create generator for dust model
gen <- odin.dust::odin_dust("models/sir_adjoint.R")

# Create dust data
incidence <- read.csv("data/incidence.csv")
incidence <- data.frame(
  time = incidence$day * 4,
  cases_observed = incidence$cases)
d <- dust::dust_data(incidence)

# Create set of parameters
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)

# Create new deterministic model with data attached
mod <- gen$new(pars, 0, 1, deterministic = TRUE)
mod$set_data(d)

compute_gradient <- function(mod, theta, trans, pd_trans, t = 0){
  pars <- as.list(trans(theta))
  mod$update_state(pars, time = t)
  res_adj <- mod$run_adjoint()
  list(log_likelihood = res_adj$log_likelihood,
       gradient = pd_trans(theta)*res_adj$gradient)
}

g <- function(theta) {exp(theta)}
dg <- function(theta) {exp(theta)}

HMC_step <- function(mod, current_theta, epsilon, L, g, dg){
  #browser()
  current_v <- rnorm(length(theta),0,1) # independent standard normal variates
  theta <- current_theta
  v <- current_v

  # Make a half step for momentum at the beginning
  v <- v - epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2

  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    theta <- theta + epsilon * v
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) v <- v - epsilon * compute_gradient(mod, theta, g, dg)$gradient }

  # Make a half step for momentum at the end.
  v <- v - epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  v <- -v
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- compute_gradient(mod, current_theta, g, dg)$log_likelihood
  current_K <- sum(current_v^2) / 2
  proposed_U <- compute_gradient(mod, theta, g, dg)$log_likelihood
  proposed_K <- sum(v^2) / 2
  if(is.na(exp(current_U-proposed_U+current_K-proposed_K))) { browser()
    }
  # Accept or reject the state at end of trajectory, returning either # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (theta) # accept
  } else {
    return (current_theta) # reject
  }
}

theta <- log(unlist(pars))
n_steps <- 1000
theta_chain <- NULL
for(i in seq(n_steps)){
  theta <- HMC_step(mod, theta, 0.0001, 10, g, dg)
  theta_chain <- rbind(theta_chain, theta)
}

plot(theta_chain[,"gamma"], type="l")


