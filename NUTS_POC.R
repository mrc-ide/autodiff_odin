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
  pars <- trans(theta)
  mod$update_state(pars, time = t)
  res_adj <- mod$run_adjoint()
  list(log_likelihood = -res_adj$log_likelihood,
       gradient = -pd_trans(theta)*res_adj$gradient)
}

hamiltonian <- function(theta, r, mod, g, dg){
  #Note that we only need the likelihood here so no need to also calculate
  #the gradient as we do here
  sum(r^2) / 2  - compute_gradient(mod, theta, g, dg)$log_likelihood
}

# perform 1 leafrog integration of step epsilon
leapfrog <- function(mod, current_theta, current_r, epsilon, g, dg){
  # initialise to the current value of theta and r
  theta <- current_theta
  r <- current_r
  # Make a half step for momentum
  r <- r - epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2
  # Make a full step for theta
  theta <- theta + epsilon * r
  # Make a half step for momentum
  r <- r - epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2
  return(list(theta = theta, r = r))
}

# function from the NUTS paper
# the NUTS paper fix the initial epsilon value to 1,
# in practice this can lead to very big leaps
# here it is given by user
# and generate NaN
find_epsilon1 <- function(mod, theta, g, dg, init_eps){
  epsilon <- init_eps
  current_theta <- theta
  current_r <- rnorm(length(theta),0,1)
  theta_r_prop <- leapfrog(mod, theta, current_r, epsilon, g, dg)
  if(exp(hamiltonian(theta_r_prop$theta,
                     theta_r_prop$r, mod, g, dg) -hamiltonian(current_theta,
                                                              current_r, mod, g, dg)) > 0.5) a <- 1 else a <- -1
  while(exp(a*(hamiltonian(theta_r_prop$theta,
                           theta_r_prop$r, mod, g, dg) -hamiltonian(current_theta,
                                                                    current_r, mod, g, dg))) > 2^-a)
  {
    epsilon <- 2^a*epsilon
    theta_r_prop <- leapfrog(mod, theta, current_r, epsilon, g, dg)
  }
  epsilon
}

build_tree <- function(theta, r, u, v, j, epsilon, theta_0, r_0, mod, g, dg, delta = 1000)
{
  if(j==0) {
    # base case, one leapfrog in direction v
    theta_r_prop <- leapfrog(mod, theta, r, epsilon, g, dg)
    H_prop <- hamiltonian(theta_r_prop$theta, theta_r_prop$r, mod, g, dg)
    H_0 <- hamiltonian(theta_0, r_0, mod, g, dg)
    n <- u <= exp(H_prop)
    s <- u < exp(delta + H_prop)
    return(list(theta_minus = theta_r_prop$theta,
           r_minus = theta_r_prop$r,
           theta_plus = theta_r_prop$theta,
           r_plus = theta_r_prop$r,
           theta_prop = theta_r_prop$theta,
           n_prop = n,
           s_prop = s,
           alpha = min(1, exp(H_prop-H_0)),
           n_alpha = 1)
    )
  }
}

g <- function(theta) {as.list(exp(theta))}
dg <- function(theta) {exp(theta)}
theta <- log(unlist(pars))
epsilon <- find_epsilon1(mod, theta, g, dg, 0.0001)
r <- rnorm(length(theta),0,1)
u <- runif(1)*exp(hamiltonian(theta, r, mod, g, dg))
tree <- build_tree(theta, r, u, v=-1, j=0, epsilon, theta, r, mod, g, dg, delta = 1000)
