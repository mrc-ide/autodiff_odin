#remotes::install_github("mrc-ide/odin@mrc-4277", upgrade = FALSE, force = TRUE)
#remotes::install_github("mrc-ide/odin.dust@mrc-4278", upgrade = FALSE, force = TRUE)
#remotes::install_github("mrc-ide/dust", upgrade = FALSE, force = TRUE)

generator <- odin.dust::odin_dust("models/logistic_growth_normal_obs.R")

mod <- generator$new(pars=list(r=1, sd_noise=5), time = 0, n_particles = 1)

n_obs <- 20
tt <- seq(0, 25, length.out = 101)
t_obs <- seq(1, n_obs)

y <- mod$simulate(tt)[,1,]
plot(tt, y, xlab = "Time", ylab = "N", main = "", type="l", ylim= c(0,120))

#default parameters
N0 <- 1
K <- 100
r <- 1

#computes the mid-point time (t0 such that N(t0)=K/2)
t0 <- log(K/N0-1)/r

#plot the analytical solution of the black model
N_obs <- K/(1+exp(-r*(t_obs-t0)))
points(t_obs, N_obs)

#generates noisy observations/data
sd_noise <- 5

d_df <- data.frame(time = t_obs,
                   observed = rnorm(N_obs,N_obs, sd_noise))

d <- dust::dust_data(d_df)

points(t_obs, d_df$observed, pch=19, col="grey")

#reversing the ode's
generator_reverse <- odin.dust::odin_dust("models/reverse_AD_logistic.R")

reverse_mod <- generator_reverse$new(pars= list(r=1, N_end=N_obs[20], t_end=20), time=0, n_particles = 1)

tt <- seq(0, 25, length.out = 101)
reverse_y <- reverse_mod$simulate(tt)[,1,]
lines(reverse_y[2,], reverse_y[1,], lty=3, col="green")

##As we can see this is slightly unstable!!!
##plot the error
plot(reverse_y[2,81:1], reverse_y[1,81:1]- y[1:81])

#Let's try to do the same now by refreshing the state at each observation point
for(i in seq_along(t_obs[-1])){
  obs_number <- length(t_obs)-i+1
  reverse_mod$set_user(r=1, N_end=N_obs[obs_number], t_end=t_obs[obs_number])
  reverse_between_obs <- reverse_mod$run(tt[tt<=t_obs[obs_number]-t_obs[obs_number-1]])
}

#Gradient of the analytical solution
ff <- parse(text="K/(1+exp(-r*t)*(K/N0-1))")
D(ff, "K")
D(ff, "N0")
D(ff, "r")

ff <- parse(text="r*N*(1-N/K)")
ff2 <- parse(text="r*N-r*N^2/K")
D(ff, "N")
D(ff2, "N")
D(ff, "K")
D(ff, "r")
D(ff, "N0")

