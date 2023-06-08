#odin SIR model
sir <- odin.dust::odin_dust("models/sir_4_AD.R")
beta <- 0.25
gamma <- 0.1
I0 <- 1
pars <- list(beta = beta, gamma = gamma, I0 = I0)
sir_model <- sir$new(pars, 0, 1)
y <- sir_model$simulate(seq(0,400, by=4))

#Plot of the model output (daily incidence)
plot(y[6,1,-1])

#Compiling the stan model
m <- rstan::stan_model("models/SIR.stan")

incidence <- read.csv("data/incidence.csv")

#Building the stanfit object
stan_fit <- rstan::stan(file = "models/SIR.stan", data = list(
  T=100,
  Y=incidence$cases,
  freq=4), chains = 0)

#Running the gradient
rstan::grad_log_prob(stan_fit,
                     rstan::unconstrain_pars(stan_fit,list(beta=0.25, gamma = 0.1, I0=1)))

#Checking that stan output the same than the odin model
sim <- rstan::constrain_pars(stan_fit,
                      rstan::unconstrain_pars(stan_fit,list(beta=0.25, gamma = 0.1, I0=1)))

lines(sim$Y_model, col="red")
