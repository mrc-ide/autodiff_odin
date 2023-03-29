incidence <- read.csv("data/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

sir <- odin.dust::odin_dust("models/sir_4_AD.R")
adj_sir <- odin.dust::odin_dust("models/adjoint_sir_4_AD.R")

compare <- function(state, observed, pars = NULL) {
  modelled <- state["incidence", , drop = TRUE]
  lambda <- modelled #+ rexp(length(modelled), 1e6)
  dpois(observed$cases, lambda, log = TRUE)
}

index <- function(info) {
  list(run = c(incidence = info$index$cases_inc),
       state = c(t = info$index$time,
                 I = info$index$I,
                 cases = info$index$cases_inc))
}

data_input <- NULL
for(i in 1:100){
  data_input <- c(data_input,rep(0,3),data$cases[i])
}

num_grad_sir <- function(beta = 0.25, gamma = 0.1, filter, h = 1e-6){
  LL_c <- filter$run(list(beta = beta, gamma = gamma))
  LL_beta_h <- filter$run(list(beta = beta+h, gamma = gamma))
  LL_gamma_h <- filter$run(list(beta = beta, gamma = gamma+h))
  (c(LL_beta_h,LL_gamma_h)-LL_c)/h
}

beta <- 0.25
gamma <- 0.1

pars <- list(beta = beta, gamma = gamma, I0 = 1)
sir_model <- sir$new(pars, 0, 1)
y <- sir_model$simulate(seq(0,400))
pars <- list(beta = beta,
             gamma = gamma,
             I0 = 1,
             main_states = y[,1,],
             data_input = data_input,
             total_steps = 400)
adj_model <- adj_sir$new(pars, 0, 1)
adj_y <- adj_model$simulate(seq(0,400))
adj_y[7:8,1,401]

#Create the filter
filter <- mcstate::particle_filter$new(data, model = sir, n_particles = 1,
                                       compare = compare, index = index)

#Run it
f_res <- filter$run(list(beta = beta, gamma = gamma, I0 = 1), save_history = TRUE)
xx_h <- filter$history()[3,1,-1]

#Get the derivative through numerical differentiation
num_grad_sir(beta,gamma,filter)

#Check it using the data_input vector and model states
kk<-data_input[rep(c(rep(FALSE,3),TRUE),100)]
xx <- y[6,1,c(FALSE,rep(c(rep(FALSE,3),TRUE),100))]
sum(-xx+kk*log(xx)-log(factorial(kk)))
sum(-xx_h+kk*log(xx_h)-log(factorial(kk)))

#Compare dll/dD
-1+kk/xx
