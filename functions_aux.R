#Convert Cholesky lower triangular matrix into VCV matrix
Chol2Cov <- function(x){
  x %*% t(x)
}

#Compute the log-posterior of the parameter using the filter and parameter characteristics provided
#Note that theta is in [-inf,inf]^2 while x in the "epidemiological space" is in [0,+inf]^2
logP <-function(theta, mcmc_pars, filter){
  x <- exp(theta)
  LL <- filter$run(save_history = TRUE, pars = list(dt = dt, beta = x[1], gamma = x[2]))
  LPr <- mcmc_pars$prior(x)

  LL+LPr
}

#Computes posterior for several samples
calculate_posterior_map <- function(parameter_samples, mcmc_pars, filter){
  apply(parameter_samples, 1, function(x)
    logP(x,mcmc_pars, filter))
}

#function to compute the gradient of the Log posterior
gradient_LP <- function(theta, mcmc_pars, filter, eps = 1e-4){
  n <- length(theta)
  theta_h <- diag(eps,n) + matrix(rep(theta,n),ncol = n)

  LP_theta <- logP(theta, mcmc_pars, filter)

  LL_h <- calculate_posterior_map(t(theta_h), mcmc_pars, filter)

  list(LP = LP_theta,
       grad_LP = (LL_h-LP_theta)/eps)
}

gradient_sir_with_LL <- function(x, sir_gen, adj_gen, data_input){
  beta <- x[1]
  gamma <- x[2]
  I0 <- 10 #x[3]
  pars <- list(beta = beta, gamma = gamma, I0 = I0)
  sir_model <- sir_gen$new(pars, 0, 1)
  y <- sir_model$simulate(seq(0,400))
  pars <- list(beta = beta,
               gamma = gamma,
               I0 = I0,
               main_states = y[,1,],
               data_input = data_input,
               total_steps = 400)
  adj_model <- adj_gen$new(pars, 0, 1)
  adj_y <- adj_model$simulate(400)
  list(LL = filter$run(pars = list(dt = 0.25,
                                    beta = beta,
                                    gamma = gamma)),
       grad_LL = adj_y[7:9,1,1])
}
