# remotes::install_github("mrc-ide/odin@mrc-4358", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/odin.dust@mrc-4359", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/dust@mrc-4307", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/mcstate", upgrade = FALSE, force = TRUE)

# Create generator for dust model
gen <- odin.dust::odin_dust("models/sir_adjoint.R")

# Create dataset
incidence <- read.csv("data/incidence.csv")
incidence <- data.frame(
  time = incidence$day * 4,
  cases_observed = incidence$cases)

# Create set of parameters
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)

# Set up the data for the particle filter
d_df <- incidence
d_df$t <- 1:100 #rename to respect convention that time should not be called time ;)
pf_data <- mcstate::particle_filter_data(d_df, "t", rate=4, initial_time = 0)

# Creating the filter
filter <- mcstate::particle_deterministic$new(data=pf_data, gen, compare = NULL)

# Running the filter for likelihood estimation
filter$run(pars = pars)

# Run mcmc
beta <- mcstate::pmcmc_parameter("beta", 0.2, min = 0)
gamma <- mcstate::pmcmc_parameter("gamma", 0.1, min = 0, prior = function(p)
  dgamma(p, shape = 1, scale = 0.2, log = TRUE))
I0 <- mcstate::pmcmc_parameter("I0", 1, min = 0)

proposal_matrix <- diag(0.1, 2)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta = beta, gamma = gamma),
                                           proposal_matrix)
n_steps <- 500000
control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE)
pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)
plot(log(pmcmc_run$pars[,1]),log(pmcmc_run$pars[,2]))

# Set up the data to attach to dust model
d <- dust::dust_data(incidence)

# Create new deterministic model with data attached
mod <- gen$new(pars, 0, 1, deterministic = TRUE)
mod$set_data(d)

# Running the likelihood evaluation using the odin model, compare and data
mod$run_adjoint()

# Simulating the model and calculating the likelihood based on this
mod$update_state(time = 0)
y <- mod$simulate(c(1:100)*4)
sum(dpois(x = incidence$cases_observed, lambda = y[mod$info()$index$cases_inc,1,], log = TRUE))

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
  print(current_theta)
  print(current_r)
  print(epsilon)
  theta <- current_theta
  r <- current_r
  # Make a half step for momentum
  r <- r + epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2
  # Make a full step for theta
  theta <- theta + epsilon * r
  # Make a half step for momentum
  r <- r + epsilon * compute_gradient(mod, theta, g, dg)$gradient / 2
  return(list(theta = theta, r = r))
}

# function from the NUTS paper
# the NUTS paper fix the initial epsilon value to 1,
# in practice this can lead to very big leaps
# here it is given by user
# and generate NaN
find_epsilon1 <- function(mod, theta, g, dg, init_eps){
  #browser()
  epsilon <- init_eps
  r <- rnorm(length(theta),0,1)
  theta_r_prop <- leapfrog(mod, theta, r, epsilon, g, dg)
  if(exp(hamiltonian(theta, r, mod, g, dg)-
         hamiltonian(theta_r_prop$theta,
                     theta_r_prop$r, mod, g, dg)) > 0.5) a <- 1 else a <- -1
  print(paste0("a= ", a))
  while(exp(a*(hamiltonian(theta,
                           r, mod, g, dg) -hamiltonian(theta_r_prop$theta,
                             theta_r_prop$r, mod, g, dg))) > 2^-a)
  {
    print(paste0(epsilon/init_eps, "-> alpha: ",
                 exp(a*(hamiltonian(theta, r, mod, g, dg) -
                          hamiltonian(theta_r_prop$theta,
                                      theta_r_prop$r, mod, g, dg)))))
    epsilon <- 2^a*epsilon
    theta_r_prop <- leapfrog(mod, theta, r, epsilon, g, dg)
  }
  print(paste0(epsilon/init_eps, "-> alpha: ",
               exp(a*(hamiltonian(theta, r, mod, g, dg) -
                        hamiltonian(theta_r_prop$theta,
                                    theta_r_prop$r, mod, g, dg)))))
  print(paste0("eps= ",epsilon, " --- (theta,r) = [", theta_r_prop$theta, ", ", theta_r_prop$r, "]"))
  epsilon
}

build_tree <- function(theta, r, u, v, j, epsilon, theta_0, r_0, mod, g, dg, delta = 1000)
{
  #browser()
  if(j==0) {
    #browser()
    # base case, one leapfrog in direction v
    theta_r_prop <- leapfrog(mod, theta, r, v*epsilon, g, dg)
    H_prop <- hamiltonian(theta_r_prop$theta, theta_r_prop$r, mod, g, dg)
    H_0 <- hamiltonian(theta_0, r_0, mod, g, dg)
    n <- as.integer(u <= exp(H_prop))
    lines(c(theta["beta"],theta_r_prop$theta["beta"]),c(theta["gamma"],theta_r_prop$theta["gamma"]))
    points(theta_r_prop$theta["beta"],theta_r_prop$theta["gamma"])
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
  } else { #j>0
    result_list <- build_tree(theta, r, u, v, j-1, epsilon, theta_0, r_0, mod, g, dg, delta)
    #browser()
    if(result_list$s_prop){ # continue the tree unless stop condition is reached
      if(v==-1){
        alternative_list <- build_tree(result_list$theta_minus, result_list$r_minus,
                               u, v, j-1, epsilon, theta_0, r_0, mod, g, dg, delta)
        result_list$theta_minus <- alternative_list$theta_minus
        result_list$r_minus <- alternative_list$r_minus
      } else { #v==1
        alternative_list <- build_tree(result_list$theta_plus, result_list$r_plus,
                               u, v, j-1, epsilon, theta_0, r_0, mod, g, dg, delta)
        result_list$theta_plus <- alternative_list$theta_plus
        result_list$r_plus <- alternative_list$r_plus
      }
      sum_n_prop <- result_list$n_prop+alternative_list$n_prop
      if(sum_n_prop > 0)
        if(runif(1)<alternative_list$n_prop/sum_n_prop)
          result_list$theta_prop <- alternative_list$theta_prop
      result_list$alpha <- result_list$alpha + alternative_list$alpha
      result_list$n_alpha <- result_list$n_alpha + alternative_list$n_alpha
      result_list$s_prop <- alternative_list$s_prop &
        ((result_list$theta_plus-result_list$theta_minus)%*%result_list$r_minus >= 0) &
        ((result_list$theta_plus-result_list$theta_minus)%*%result_list$r_plus >= 0)
      result_list$n_prop <- sum_n_prop
    }
    return(result_list)
  }
}

g <- function(theta) {as.list(exp(theta))}
dg <- function(theta) {exp(theta)}
theta <- log(unlist(pars))
plot(theta["beta"],theta["gamma"],
     xlim=c(theta["beta"]-2,theta["beta"]+2),
     ylim=c(theta["gamma"]-2,theta["gamma"]+2), pch=19, col="red")
epsilon <- find_epsilon1(mod, theta, g, dg, 0.000001)
for(i in 1:100){
  r <- rnorm(length(theta),0,1)
  u <- runif(1)*exp(hamiltonian(theta, r, mod, g, dg))
  tree <- build_tree(theta, r, u, v=1, j=12, epsilon/10, theta, r, mod, g, dg, delta = 1000)
  print(tree$n_prop)
}

theta <- c(-106.30584, 100.80329, -21.51924)
r <- rnorm(length(theta),0,1)
leapfrog(mod, theta, r, epsilon, g, dg)

# plot(theta["beta"],theta["gamma"],
#      xlim=c(theta["beta"]-3,theta["beta"]+2),
#      ylim=c(theta["gamma"]-2,theta["gamma"]+2), pch=19, col="red")
#
# theta0 <- log(unlist(pars))
# M <- 1000
# M_adapt <- 100
# D_max <- 1000
# epsilon0 <- find_epsilon1(mod, theta, g, dg, 0.0001)
# #epsilon0 <- 0.005
# mu <- log(10*epsilon0)
# theta_m <- matrix(rep(theta0, M+1), ncol = length(theta0), byrow = TRUE)
# colnames(theta_m) <- names(theta0)
#
# tt <- rbind(theta0 ,theta0 ,theta0 )
#
# for(i in 1:M)
# {
#   r0 <- rnorm(length(theta),0,1)
#   print(r0)
#   u <- runif(1)*exp(hamiltonian(theta_m[i,], r0, mod, g, dg))
#   theta_minus <- theta_m[i,]
#   theta_plus <- theta_m[i,]
#   r_minus <- r0
#   r_plus <- r0
#   j <- 0
#   theta_prop <- theta_m[i,]
#   n <- 1
#   s <- TRUE
#   while(s){
#     v <- sample(c(-1,1),1)
#     if(v==-1){
#       tree_list <- build_tree(theta_minus, r_minus,
#                                      u, v, j, epsilon0, theta_m[i,], r0, mod, g, dg, D_max)
#     } else {
#       tree_list <- build_tree(theta_plus, r_plus,
#                               u, v, j, epsilon0, theta_m[i,], r0, mod, g, dg, D_max)
#     }
#     #browser()
#     if(tree_list$s_prop)
#       if(runif(1)<min(1,tree_list$n_prop/n))
#         theta_m[i+1,] <- tree_list$theta_prop
#     n <- n + tree_list$n_prop
#     s <- tree_list$s_prop &
#       ((tree_list$theta_plus-tree_list$theta_minus)%*%tree_list$r_minus >= 0) &
#       ((tree_list$theta_plus-tree_list$theta_minus)%*%tree_list$r_plus >= 0)
#     j <- j+1
#   }
# }



