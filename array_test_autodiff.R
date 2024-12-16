# odin model
sir <- odin2::odin({

  # Equations for transitions between compartments by age group
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  update(incidence) <- incidence + n_SI[1] + n_SI[2]

  # Individual probabilities of transition:

  p_SI[] <- 1 - exp(-lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R

  # Calculate force of infection

  # age-structured contact matrix: m[i, j] is mean number of contacts an
  # individual in group i has with an individual in group j per time unit

  m <- parameter()

  # here s_ij[i, j] gives the mean number of contacts and individual in group
  # i will have with the currently infectious individuals of group j
  s_ij[, ] <- m[i, j] * I[j]

  # lambda[i] is the total force of infection on an individual in group i
  lambda[] <- beta * (s_ij[i, 1] + s_ij[i, 2])

  # Draws from binomial distributions for numbers changing between
  # compartments:

  n_SI[] <- S[i] * p_SI[i]
  n_IR[] <- I[i] * p_IR

  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(incidence, zero_every = 1) <- 0

  # User defined parameters - default in parentheses:
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.000165)
  gamma <- parameter(0.1)
  rho <- parameter(0.1)

  # Dimensions of arrays
  dim(S0) <- 2
  dim(I0) <- 2
  dim(S) <- 2
  dim(I) <- 2
  dim(R) <- 2
  dim(n_SI) <- 2
  dim(n_IR) <- 2
  dim(p_SI) <- 2
  dim(m) <- c(2, 2)
  dim(s_ij) <- c(2, 2)
  dim(lambda) <- 2

  cases <- data()
  cases ~ Poisson(rho * incidence)
})

# Generating data
pars <- list(
  S0 = c(1000,1000),
  I0 = c(10,0),
  m = matrix(c(1,0.5,1,0.5), ncol = 2))
sys <- dust2::dust_system_create(sir, pars, dt = 0.25)
dust2::dust_system_set_state_initial(sys)
dust2::dust_system_state(sys)
t <- seq(1, 100)
y <- dust2::dust_system_simulate(sys, t)
set.seed(42)
cases <- rpois(length(t), y[7,seq_along(t)] * 0.1)
data <- data.frame(time = t, cases = cases)

# Filter and likelihood
filter <- dust2::dust_unfilter_create(sir, 0, data)

dust2::dust_likelihood_run(filter, pars)
