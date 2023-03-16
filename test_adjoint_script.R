#testing the adjoint model on 1 step and LL depending on the final state (step=1)

sir <- odin.dust::odin_dust("models/sir_4_AD.R")

#d <- data.frame(t = c(1,2), y = c(sample(10,1),NA))
d <- data.frame(t = c(1,2), y = sample(2,2))
data <- mcstate::particle_filter_data(d, "t", 1, initial_time = 0)
data <- data[-1,]

compare <- function(state, observed, pars = NULL) {
  modelled <- state["incidence", , drop = TRUE]
  lambda <- modelled #+ rexp(length(modelled), 1e6)
  dpois(observed$y, lambda, log = TRUE)
}

index <- function(info) {
  list(run = c(incidence = info$index$cases_inc),
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))
}

filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
                                       compare = compare, index = index)

pars <- list(beta = 0.25, gamma = 0.1, I0 = 1, freq = 1)

LL1 <- filter$run(pars)

#testing the pipeline with running the model followed by calculating the poisson density

mod <- sir$new(pars, 0, 1)
y <- mod$simulate(c(0,1))
LL2 <- dpois(data$y[1],y[5,1,2], log = TRUE)

LL_beta <- function(beta){
  pars <- list(beta = beta, gamma = 0.1, I0 = 1, freq = 1)
  mod$initialize(pars, 0, 1)
  y <- mod$simulate(c(0,1))
  dpois(data$y[1],y[5,1,2], log = TRUE)
}

#plot the LL surface as a function of beta for gamma=0.1 and I0=1
beta_seq <- seq(from=.1, to=10, length.out=100)
res <- sapply(beta_seq, LL_beta)
plot(beta_seq, res, type="l")

mod$info()$index

#Calculation of the intermediate variable to go from y[,1,1] to y[,1,2]
gamma <- 0.1
beta <- 0.25

S <- y[2,1,1]
R <- y[3,1,1]
I <- y[4,1,1]
cases_cumul <- y[5,1,1]
cases_inc <- y[6,1,1]

p_IR <- 1 - exp(-(gamma))
dt <- 0.25

N <- S + I + R
p_inf <- beta * I / N
p_SI <- 1 - exp(-(p_inf))
n_SI <- S * p_SI * dt
n_IR <- I * p_IR * dt

S_plus <- S - n_SI
I_plus <- I + n_SI - n_IR
R_plus <- R + n_IR
cases_cumul_plus <- cases_cumul + n_SI
cases_inc_plus <- if (step %% freq == 0) n_SI else cases_inc + n_SI

LL3 <- dpois(data$y[1],cases_cumul_plus, log = TRUE)

#All three methods give the same result for the pois LL

#Let's try to get the value dLL/dbeta at beta = 0.25, gamma = 0.1
der1 <- (LL_beta(0.25+1e-6)-LL_beta(0.25))/1e-6

#plotting the tangent
abline(LL3 - der1*0.25, der1, col="red")

sir <- odin.dust::odin_dust("models/sir_4_AD.R")
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1, freq = 1)
mod <- sir$new(pars, 0, 1)
y <- mod$simulate(c(0,1))
m_y <- y[,1,]
data_input <- c(-1, data$y[1]) #rep(data$y[1],2)
adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")
pars_adj <- list(beta = 0.25,
             gamma = 0.1,
             I0 = 1,
             freq = 1,
             main_states = m_y,
             data_input = data_input,
             total_steps = 1)
adj_mod <- adj_sir$new(pars_adj, 0, 1)
adj_y <- adj_mod$simulate(c(0,1))

#Calculate the LL and the perturbated LL
LL_m <- LL_beta(0.25)
LL_p <- LL_beta(0.25+1e-6)
LL_p-LL_m

gamma <- 0.1
beta_m <- 0.25
beta_p <- 0.25+1e-6

S <- y[2,1,1]
R <- y[3,1,1]
I <- y[4,1,1]
cases_cumul <- y[5,1,1]
cases_inc <- y[6,1,1]

p_IR <- 1 - exp(-(gamma))
dt <- 1

N <- S + I + R
p_inf_m <- beta_m * I / N
p_inf_p <- beta_p * I / N
p_SI_m <- 1 - exp(-(p_inf_m))
p_SI_p <- 1 - exp(-(p_inf_p))
n_SI_m <- S * p_SI_m * dt
n_SI_p <- S * p_SI_p * dt
n_IR <- I * p_IR * dt

#time_plus <- (step + 1) * dt
S_plus_m <- S - n_SI_m
S_plus_p <- S - n_SI_p
I_plus_m <- I + n_SI_m - n_IR
I_plus_p <- I + n_SI_p - n_IR
R_plus <- R + n_IR
cases_inc_plus_m <- cases_inc + n_SI_m
cases_inc_plus_p <- cases_inc + n_SI_p
#cases_inc_plus <- if (step %% freq == 0) n_SI else cases_inc + n_SI


