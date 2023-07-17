# remotes::install_github("mrc-ide/odin", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/dust", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/mcstate", upgrade = FALSE, force = TRUE)

generator <- odin.dust::odin_dust("models/logistic_growth_normal_obs.R")

#default parameters
N0 <- 2
K <- 100
r <- 1

mod <- generator$new(pars=list(r=r, N0=N0, K=K, sd_noise=5), time = 0, n_particles = 1)

n_obs <- 20
t_obs <- seq(1, n_obs)

#Plot the trajectory
tt <- seq(0, 25, length.out = 101)
y <- mod$simulate(tt)[,1,]
plot(tt, y, xlab = "Time", ylab = "N", main = "", type="l", ylim= c(0,120))

#computes the mid-point time (t0 such that N(t0)=K/2)
t0 <- log(K/N0-1)/r

#plot the analytical solution of the plotted model
N_obs <- K/(1+exp(-r*(t_obs-t0)))
points(t_obs, N_obs)

#generates noisy observations/data
sd_noise <- 5
d_df <- data.frame(time = t_obs,
                   observed = rnorm(N_obs,N_obs, sd_noise))
points(t_obs, d_df$observed, pch=19, col="grey")


#Example of calculate the contribution to the log-likelihood of one data point
#using the inline compare and the generated data
d <- dust::dust_data(d_df)
mod$initialize(pars=list(r=1, sd_noise=5), time = 0, n_particles = 1)
mod$set_data(d)
yy <- mod$run(20)
mod$compare_data()

#Run system until last observation
last_obs <- length(t_obs)
mod$initialize(pars=list(r=1, sd_noise=5), time = 0, n_particles = 1)
mod$set_data(d)
y_end <- mod$run(last_obs)

#reversing the ode's
generator_reverse <- odin.dust::odin_dust("models/reverse_AD_logistic.R")

N_curr <- N_obs[last_obs]
t_curr <- t_obs[last_obs]
adj_N_curr <- 0
adj_K_curr <- 0
adj_r_curr <- 0
y_curr <- y_end

reverse_mod <- generator_reverse$new(pars= list(r=1,
                                                N_end=N_curr,
                                                t_end=t_curr,
                                                adj_N_end=adj_N_curr,
                                                adj_K_end=0,
                                                adj_r_end=0), time=0, n_particles = 1)

contribution_data <- function(observed, model, sd){
  (observed-model)/sd^2
}

for(i in seq_along(t_obs)[-n_obs]){
  adj_N_curr <- adj_N_curr + contribution_data(N_curr, y_curr, sd_noise)
  reverse_mod$initialize(pars= list(r=1,
                                    N_end=N_curr,
                                    t_end=t_curr,
                                    adj_N_end=adj_N_curr,
                                    adj_K_end=adj_K_curr,
                                    adj_r_end=adj_r_curr), time=0, n_particles = 1)
  reverse_mod$run(t_curr-t_obs[n_obs-i])

  N_curr <- reverse_mod$info()$index$N
  t_curr <- t_obs[n_obs-i]
  adj_N_curr <- reverse_mod$info()$index$adj_N
  adj_K_curr <- reverse_mod$info()$index$adj_K
  adj_r_curr <- reverse_mod$info()$index$adj_r
}

reverse_mod$info()$index

##As we can see this is slightly unstable!!!
##plot the error
#plot(reverse_y[2,81:1], reverse_y[1,81:1]- y[1:81])

#Let's try to do the same now by refreshing the state at each observation point
# for(i in seq_along(t_obs[-1])){
#   obs_number <- length(t_obs)-i+1
#   reverse_mod$set_user(r=1, N_end=N_obs[obs_number], t_end=t_obs[obs_number])
#   reverse_between_obs <- reverse_mod$run(tt[tt<=t_obs[obs_number]-t_obs[obs_number-1]])
# }

#Gradient of the analytical solution
ff <- parse(text="K/(1+exp(-r*t)*(K/N0-1))")
D(ff, "K")
D(ff, "N0")
D(ff, "r")

ff <- parse(text="r*N*(1-N/K)")
ff2 <- parse(text="r*N-r*N^2/K")
ff3 <- parse(text="-log(sigma)-a/sigma^2")
D(ff3,"sigma")
D(ff, "N")
D(ff2, "N")
D(ff, "K")
D(ff, "r")
D(ff, "N0")

ff3 <- parse(text="-log(sigma)-1/2*((N_obs-N)/sigma)^2")
D(ff3, "N")

#Not working analysis - self-contained code

#odin_dust model
generator <- odin.dust::odin_dust("
deriv(N) <- r * N * (1 - N / K)
initial(N) <- N0

N0 <- user(1)
K <- user(100)
r <- user()

sd_noise <- user()

observed <- data()
compare(observed) ~ normal(N, sd_noise)")

#default parameters
N0 <- 2
K <- 100
r <- 1

#generating "data"
mod <- generator$new(pars=list(r=r, N0=N0, K=K, sd_noise=5), time = 0, n_particles = 1)

n_obs <- 20
t_obs <- seq(1, n_obs)
t0 <- log(K/N0-1)/r #computes the mid-point time (t0 such that N(t0)=K/2)
N_obs <- K/(1+exp(-r*(t_obs-t0))) #derive the analytical solution model
sd_noise <- 5
d_df <- data.frame(time = t_obs,
                   observed = rnorm(N_obs,N_obs, sd_noise))

#creates data for pf
d_df$t <- d_df$time #rename to respect convention
pf_data <- mcstate::particle_filter_data(d_df, "t", rate=NULL, initial_time = 0)

#cannot create the pf
filter <- mcstate::particle_filter$new(data=pf_data, generator, n_particles = 1, compare = NULL)
h <- 1e-6
(filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=5+h))-filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=5)))/h

filter$run(pars=list(r=r, N0=N0, K=K, sd_noise=5))
sum(dnorm(x = d_df$observed, mean = N_obs, sd = sd_noise, log = TRUE))

# Gradient function with respect to K
partial_derivative_K <- function(K, N0, r, t) {
  gradient <- 1 / (1 + ((K - N0) / N0) * exp(-r * t))^2
  return(gradient)
}

# Gradient function with respect to N0
partial_derivative_N0 <- function(K, N0, r, t) {
  gradient <- (-K / N0^2) * exp(-r * t) / (1 + ((K - N0) / N0) * exp(-r * t))^2
  return(gradient)
}

# Gradient function with respect to r
partial_derivative_r <- function(K, N0, r, t) {
  gradient <- (K * t * (K - N0) * exp(-r * t)) / (N0 * (1 + ((K - N0) / N0) * exp(-r * t))^2)
  return(gradient)
}

likelihood_gradient <- function(K, N0, r, sigma, N_obs, t_obs){
  #browser()
  t0 <- log(K/N0-1)/r
  g_LL_K <- 0
  g_LL_N0 <- 0
  g_LL_r <- 0
  g_LL_sigma <- 0
  for(i in seq_along(t_obs)){
    N_t <- K/(1+exp(-r*(t_obs[i]-t0)))

    #Calculate the simpler PD wrt sigma
    g_LL_sigma <- g_LL_sigma - 1/sigma + 1/sigma^3*(N_obs[i]- N_t)^2

    #Normal pdf contribution
    dlogf_dN <- (N_obs[i]- N_t)/sigma^2
    g_LL_K <- g_LL_K + dlogf_dN * partial_derivative_K(K, N0, r, t_obs[i])
    g_LL_N0 <- g_LL_N0 + dlogf_dN * partial_derivative_N0(K, N0, r, t_obs[i])
    g_LL_r <- g_LL_r + dlogf_dN * partial_derivative_r(K, N0, r, t_obs[i])

  }
  return(list(sigma = g_LL_sigma, K = g_LL_K, N0 = g_LL_N0, r = g_LL_r))
}

likelihood_gradient(K, N0, r, sd_noise, d_df$observed, t_obs)

#Run system until last observation
last_obs <- length(t_obs)
mod$initialize(pars=list(r=1, sd_noise=5), time = 0, n_particles = 1)
d <- dust::dust_data(d_df)
mod$set_data(d)

#reversing the ode's
generator_reverse <- odin.dust::odin_dust("models/reverse_AD_logistic.R")

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
                                                adj_r_end=0), time=0, n_particles = 1)

contribution_data <- function(observed, model, sd){
  (observed-model)/sd^2
}

contribution_data_sigma <- function(observed, model, sd){
  -1/sd + (observed-model)^2/sd^3
}

for(i in seq_along(t_obs)){
  y_curr <- d_df$observed[n_obs-i+1]
  adj_N_curr <- adj_N_curr + contribution_data(N_curr, y_curr, sd_noise)
  adj_sigma_curr <- adj_sigma_curr + contribution_data_sigma(N_curr, y_curr, sd_noise)
  reverse_mod$initialize(pars= list(r=1,
                                    N_end=N_curr,
                                    t_end=t_curr,
                                    adj_N_end=adj_N_curr,
                                    adj_K_end=adj_K_curr,
                                    adj_r_end=adj_r_curr), time=0, n_particles = 1)
  reverse_mod$run(t_curr-t_obs[n_obs-i])

  N_curr <- reverse_mod$state()[reverse_mod$info()$index$N]
  t_curr <- t_obs[n_obs-i]
  adj_N_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_N]
  adj_K_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_K]
  adj_r_curr <- reverse_mod$state()[reverse_mod$info()$index$adj_r]
}

reverse_mod$state()[reverse_mod$info()$index$adj_N]
reverse_mod$state()[reverse_mod$info()$index$adj_K]
reverse_mod$state()[reverse_mod$info()$index$adj_r]
adj_sigma_curr


