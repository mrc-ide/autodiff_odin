source("functions_aux.R")

sir <- odin.dust::odin_dust("models/sir_4_AD.R")

case_compare <- function(state, observed, pars = NULL) {
  incidence_modelled <- state[5, , drop = TRUE]
  incidence_observed <- observed$cases
  dpois(x = incidence_observed, lambda = incidence_modelled, log = TRUE)
}

incidence <- read.csv("data/incidence.csv")

dt <- 0.25
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt,
                                          initial_time = 0)

filter <- mcstate::particle_deterministic$new(data = sir_data,
                                              model = sir,
                                              compare = case_compare)

prior_par <- list(
  meanlog_beta = -3,
  sdlog_beta = 1.5,
  meanlog_gamma = -2,
  sdlog_gamma = 1
)

beta <- mcstate::pmcmc_parameter("beta", 0.16, min = 0,
                                 prior = function(p)
                                   dlnorm(p, meanlog = prior_par$meanlog_beta,
                                          sdlog = prior_par$sdlog_beta, log = TRUE)
)
gamma <- mcstate::pmcmc_parameter("gamma", 0.38, min = 0,
                                  prior = function(p)
                                    dlnorm(p, meanlog = prior_par$meanlog_gamma,
                                           sdlog = prior_par$sdlog_gamma, log = TRUE)
)


proposal_matrix <- matrix(c(0.0004507713, 0.001077553,
                            0.001077553, 0.002629411), nrow = 2, ncol = 2, byrow = TRUE)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta = beta, gamma = gamma), 2 * proposal_matrix)

adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")

data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

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

#Plot the mcmc samples
points(log(as.numeric(mcmc_run$pars[,"beta"])),
       log(as.numeric(mcmc_run$pars[,"gamma"])),
       xlim=c(-7,1), ylim = c(-5,1), col="red", pch=19)

control <- mcstate::pmcmc_control(
  10000,
  save_state = FALSE,
  save_trajectories = FALSE,
  progress = TRUE)
mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

plot(log(as.numeric(mcmc_run$pars[,"beta"])), type="l")

length(unique(mcmc_run$pars[,"beta"]))/length(mcmc_run$pars[,"beta"])

#Set the initial mvnorm approximation
initial_Chol_t <- matrix(c(.6,-0.1,0,.3), nrow = 2)
initial_mu_t <- c(-5,-1)

#Learning rate
rho <- 0.0001

n_iteration <- 1000
KL_stoch <- rep(0,n_iteration)
mu_t_chain <- initial_mu_t

Chol_t <- initial_Chol_t
mu_t <- initial_mu_t
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
}

posterior_value_samples <- calculate_posterior_map(draw_prior,
                                                   mcmc_pars,
                                                   filter)
n_pal <- 20
rbPal <- colorRampPalette(c('yellow','blue'))
plot(draw_prior, pch = 19, col=grey(.8),
     xlab="log(beta)", ylab="log(gamma)",
     main ="Variational inference")
legend("bottom",title="Ventile of LogPosterior",
       legend=c(1:n_pal),col =rbPal(n_pal),pch=20, horiz = TRUE, bty = "n", cex = .6)

#plot points with the colours
finite_prior <- posterior_value_samples!=-Inf
colZ <- rbPal(n_pal)[as.numeric(
  cut(log(abs(posterior_value_samples[finite_prior])),breaks = n_pal))]
points(draw_prior[finite_prior,],pch = 20,col = colZ)

lines(ellipse::ellipse( Chol2Cov(initial_Chol_t) , centre = initial_mu_t) , col='orange', lwd=3)
points(log(as.numeric(mcmc_run$pars[,"beta"])),log(as.numeric(mcmc_run$pars[,"gamma"])), xlim=c(-7,1), ylim = c(-5,1), col="red", pch=19)
lines(mu_t_chain, col="cyan")

lines(ellipse::ellipse( Chol2Cov(Chol_t) , centre = mu_t) , col='green', lwd=3)

gradient_LP(theta, mcmc_pars, filter)

gradient_sir(exp(theta), sir, adj_sir, data_input)



