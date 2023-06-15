#remotes::install_github("mrc-ide/odin@mrc-4277", upgrade = FALSE, force = TRUE)
#remotes::install_github("mrc-ide/odin.dust@mrc-4278", upgrade = FALSE, force = TRUE)

generator <- odin.dust::odin_dust("models/logistic_growth_normal_obs.R")

mod <- generator$new(r=1)

n_obs <- 20
tt <- seq(0, 25, length.out = 101)
t_obs <- seq(1, n_obs)

y <- mod$run(tt)
y2 <- mod$run(tt, 50)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "", ylim= c(0,120))
lines(y2, col = "red")

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



