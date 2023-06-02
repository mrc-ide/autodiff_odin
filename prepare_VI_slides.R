source("functions_aux.R")

######################################################
#
#      Set up the data, model and filter
#
######################################################

#creates the generator for the SIR dust model
sir <- odin.dust::odin_dust("models/sir_4_AD.R")

#load the incidence data for the simple SIR example
incidence <- read.csv("data/incidence.csv")

#set the timestep
dt <- 0.25

#put the data in the format needed for the filter
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt,
                                          initial_time = 0)

#creates the compare function
case_compare <- function(state, observed, pars = NULL) {
  incidence_modelled <- state[5, , drop = TRUE]
  incidence_observed <- observed$cases
  dpois(x = incidence_observed, lambda = incidence_modelled, log = TRUE)
}

#creates the particle filter
filter <- mcstate::particle_deterministic$new(data = sir_data,
                                              model = sir,
                                              compare = case_compare)

#list of parameters for the priors (set to 2 logNormal distribution)
prior_par <- list(
  meanlog_beta = -3,
  sdlog_beta = 3,
  meanlog_gamma = -2,
  sdlog_gamma = 5
)

#define the first parameter beta
beta <- mcstate::pmcmc_parameter("beta", 0.16, min = 0,
                                 prior = function(p)
                                   dlnorm(p, meanlog = prior_par$meanlog_beta,
                                          sdlog = prior_par$sdlog_beta, log = TRUE)
)

#define the second parameter gamma
gamma <- mcstate::pmcmc_parameter("gamma", 0.38, min = 0,
                                  prior = function(p)
                                    dlnorm(p, meanlog = prior_par$meanlog_gamma,
                                           sdlog = prior_par$sdlog_gamma, log = TRUE)
)

#proposal matrix for the MH randow walk
proposal_matrix <- matrix(c(0.0004507713, 0.001077553,
                            0.001077553, 0.002629411), nrow = 2, ncol = 2, byrow = TRUE)

#creates the mcmc_pars object
mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta = beta, gamma = gamma), 5 * proposal_matrix)

#control parameters for the pmcmc
control <- mcstate::pmcmc_control(
  10000,
  save_state = FALSE,
  save_trajectories = FALSE,
  progress = TRUE)

#run the (p)MCMC
#mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

#plot the chain for beta the first parameter
plot(log(as.numeric(mcmc_run$pars[,"beta"])), type="l")

#hist of the log posterior densities for the MCMC samples
hist(mcmc_run$probabilities[,"log_posterior"], breaks = 100)

#computes acceptance rate
length(unique(mcmc_run$pars[,"beta"]))/length(mcmc_run$pars[,"beta"])

######################################################
#
#      Plot of prior distribution and MCMC samples
#
######################################################

#Mean of the two distributions
mu <- c(prior_par$meanlog_beta,prior_par$meanlog_gamma)
#Build VCV matrix assuming independance between the two prior distributions
l <- matrix(c(prior_par$sdlog_beta^2,0,0,prior_par$sdlog_gamma^2), nrow=2)
#Draw sample from prior
draw_prior <- mvtnorm::rmvnorm(n = 1000, mean = mu, sigma = l)

#Plot samples from prior
plot(draw_prior, pch = 19, col=grey(.8),
     xlab="log(beta)", ylab="log(gamma)",
     main ="MCMC and prior")

#Plot 95% CI of the prior distribution
lines(ellipse::ellipse( l , centre = mu) , col='red', lwd=3, lty = 3)

posterior_value_samples <- calculate_posterior_map(draw_prior,
                                                   mcmc_pars,
                                                   filter)

#Plot the mcmc samples
points(log(as.numeric(mcmc_run$pars[,"beta"])),
       log(as.numeric(mcmc_run$pars[,"gamma"])),
       xlim=c(-7,1), ylim = c(-5,1), col="red", pch=19)

######################################################
#
#      Run Variational Inference
#
######################################################

#Set the initial mvnorm approximation
initial_Chol_t <- matrix(c(.6,-0.5,0,.3), nrow = 2)
initial_mu_t <- c(1,1)

#Learning rate
rho <- 0.0001

n_iteration <- 20
KL_stoch <- rep(0,n_iteration)
mu_t_chain <- initial_mu_t

Chol_t <- initial_Chol_t
mu_t <- initial_mu_t
produce_frame(0, draw_prior, posterior_value_samples, mcmc_pars,
              filter,
              initial_Chol_t,
              initial_mu_t,
              mcmc_run,
              Chol_t,
              mu_t,
              mu_t_chain)

for(t in 1:n_iteration){
  z <- rnorm(2)
  theta <- Chol_t %*% z + mu_t
  grad_theta <- gradient_LP(theta, mcmc_pars, filter)
  mu_t_minus <- mu_t
  mu_t <- mu_t + rho*grad_theta$grad_LP
  Chol_t <- Chol_t + rho*0.1 * (grad_theta$grad_LP %*% t(z)+diag(1/diag(Chol_t)))
  Chol_t[upper.tri(Chol_t)] <- 0
  diag(Chol_t)[diag(Chol_t) < 1e-4] <- 1e-4
  KL_stoch[t] <- grad_theta$LP + sum(log(diag(Chol_t)))
  mu_t_chain <- rbind(mu_t_chain,mu_t)

  produce_frame(t, draw_prior, posterior_value_samples, mcmc_pars,
                filter,
                initial_Chol_t,
                initial_mu_t,
                mcmc_run,
                Chol_t,
                mu_t,
                mu_t_chain)
}

png_files <- paste0("animation/VI_frame", sprintf("%06d", seq(0,n_iteration)), ".png")
av::av_encode_video(png_files, 'output.mp4', framerate = 24)
utils::browseURL('output.mp4')

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

# index <- function(info) {
#   list(run = c(incidence = info$index$cases_inc),
#        state = c(t = info$index$time,
#                  I = info$index$I,
#                  cases = info$index$cases_inc))
# }
#
# compare <- function(state, observed, pars = NULL) {
#   modelled <- state["incidence", , drop = TRUE]
#   lambda <- modelled #+ rexp(length(modelled), 1e6)
#   dpois(observed$cases, lambda, log = TRUE)
# }
#
# filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
#                                        compare = compare, index = index)

beta <- exp(theta[1])#0.25
gamma <- exp(theta[2])#0.1
I0 <- 10

# adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")
#
# data_input <- NULL
# for(i in 1:100){
#   data_input <- c(data_input,rep(0,3),incidence$cases[i])
# }

gradient_sir(c(beta,gamma,I0), sir, adj_sir, data_input)

# m <- rstan::stan_model("models/SIR.stan")
# stan_fit <- rstan::stan(file = "models/SIR.stan", data = list(
#   T=100,
#   Y=incidence$cases,
#   freq=4), chains = 0)
# rstan::grad_log_prob(stan_fit,
#                      rstan::unconstrain_pars(stan_fit,list(beta=beta, gamma = gamma, I0=I0)))


num_grad_sir(c(beta,gamma,I0), filter)

gradient_sir_with_LL(exp(theta), sir, adj_sir, data_input)

c(gradient_LP(theta, mcmc_pars, filter)$grad_LP/exp(theta))
