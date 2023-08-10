# R code for the "AutoDiff, ODEs and odin" presentation

# Libraries
library(odin.dust)

# Contributing functions
contribution_data_adjoint <- function(observed, model, sd){
  (observed - model) / sd^2
}

contribution_data_sigma <- function(observed, model, sd){
  -1/sd + (observed - model)^2 / sd^3
}

#default parameters
N0 <- 2
K <- 100
r <- 1
sd_noise <- 5

generator <- odin.dust::odin_dust("models/logistic_growth_normal_obs.R")
mod <- generator$new(pars=list(r=r, N0=N0, K=K, sd_noise=5), time = 0, n_particles = 1)

n_obs <- 1
t_obs <- seq(1, n_obs)

#computes the mid-point time (t0 such that N(t0)=K/2)
t0 <- log(K/N0-1)/r

#generates noisy observations/data
N_obs <- K/(1+exp(-r*(t_obs-t0)))
d_df <- data.frame(time = c(0,t_obs),
                   observed = c(NA,rnorm(n_obs,N_obs, sd_noise)))

last_obs <- length(t_obs)
mod$initialize(pars=list(r=r, N0=N0, K=K, sd_noise=5), time = 0, n_particles = 1)
y_end <- mod$run(t_obs[last_obs])

# Compile the reverse model
generator_reverse <- odin.dust::odin_dust("models/reverse_AD_logistic.R")

# Initialisation of the model and augmented state
N_curr <- N_obs[last_obs]
t_curr <- t_obs[last_obs]
adj_N_curr <- 0
adj_K_curr <- 0
adj_r_curr <- 0
adj_sigma_curr <- 0

reverse_mod <- generator_reverse$new(pars= list(r=1,
                                                N_end=N_curr,
                                                t_end=t_curr,
                                                adj_N_end=adj_N_curr,
                                                adj_K_end=0,
                                                adj_r_end=0),
                                     time=0, n_particles = 1)

# Loop through observations
for(i in seq_along(t_obs)){
  y_curr <- d_df$observed[n_obs - i + 2] #this is +2 rather than +1 because we have added an NA "observation" at t=0
  adj_N_curr <- adj_N_curr + contribution_data_adjoint(y_curr, N_curr, sd_noise)
  adj_sigma_curr <- adj_sigma_curr + contribution_data_sigma(y_curr, N_curr, sd_noise)

  reverse_mod$initialize(pars= list(r=1,
                                    N_end=N_curr,
                                    t_end=t_curr,
                                    adj_N_end=adj_N_curr,
                                    adj_K_end=adj_K_curr,
                                    adj_r_end=adj_r_curr),
                         time=0,
                         n_particles = 1)

  flow_time <- if(n_obs > i) t_curr - t_obs[n_obs - i] else t_curr
  reverse_mod$run(flow_time)

  N_curr <- reverse_mod$state()[reverse_mod$info()$index$N]
  t_curr <- t_obs[n_obs - i]
  adj_N_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_N]
  adj_K_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_K]
  adj_r_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_r]
}

gradient_AD <- c(
  adj_sigma_curr,
  reverse_mod$state()[reverse_mod$info()$index$adj_N],
  reverse_mod$state()[reverse_mod$info()$index$adj_K],
  reverse_mod$state()[reverse_mod$info()$index$adj_r])
name_var <- c("sigma","N0","K","r")
names(gradient_AD) <- name_var

gradient_AD

d <- dust::dust_data(d_df)
mod$set_data(d)

#creates data for pf
d_df$t <- d_df$time #rename to respect convention that time should not be called time ;)
pf_data <- mcstate::particle_filter_data(d_df, "t", rate=NULL, initial_time = 0)

#creating the filter
filter <- mcstate::particle_filter$new(data=pf_data, generator, n_particles = 1, compare = NULL)
h <- 1e-6
ND_sigma <- (filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=sd_noise+h))-filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=sd_noise)))/h
ND_N <- (filter$run(pars=list(r=r, N0=N0+h, K=K, sd_noise=sd_noise))-filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=sd_noise)))/h
ND_K <- (filter$run(pars=list(r=r, N0=N0, K=K+h, sd_noise=sd_noise))-filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=sd_noise)))/h
ND_r <- (filter$run(pars=list(r=r+h, N0=N0, K=K, sd_noise=sd_noise))-filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=sd_noise)))/h

gradient_ND <- c(ND_sigma,ND_N,ND_K,ND_r)

gradient_AD
gradient_ND

