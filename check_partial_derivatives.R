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
h <- 1e-6

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

####Check partial derivatives

num_PD <- function(y, step, line, h=1e-6){
  total_steps <- dim(y)[3]-1
  main_step <- step

  #perturbation vector
  hv <- rep(0,6)
  hv[line] <- h

  #y_h
  y_h <- y[6,1,1:main_step]

  #pull the states of the main model
  time <- y[1,1,main_step+1] + hv[1]
  S <- y[2,1,main_step+1] + hv[2]
  R <- y[3,1,main_step+1] + hv[3]
  I <- y[4,1,main_step+1] + hv[4]
  cases_cumul <- y[5,1,main_step+1] + hv[5]
  cases_inc <- y[6,1,main_step+1] + hv[6]

  #browser()

  #update y_h
  y_h <- c(y_h,cases_inc)

  for(i in (main_step+1):total_steps){
    #Recalculate the intermediate variables between the state[main_step] and state[main_step+1]
    p_IR <- 1 - exp(-(gamma))
    freq <- 4
    dt <- 1.0 / freq

    N <- S + I + R
    p_inf <- beta * I / N
    p_SI <- 1 - exp(-(p_inf))
    n_SI <- S * p_SI * dt
    n_IR <- I * p_IR * dt

    step <- total_steps - i - 1
    time <- (step + 1) * dt
    #print(time)
    S <- S - n_SI
    I <- I + n_SI - n_IR
    R <- R + n_IR
    cases_cumul <- cases_cumul + n_SI
    cases_inc <- if ((i-1) %% freq == 0) n_SI else cases_inc + n_SI
    y_h <- c(y_h,cases_inc)
  }

  #browser()

  kk<-data_input[rep(c(rep(FALSE,3),TRUE),100)]
  xx_h <- y_h[c(FALSE,rep(c(rep(FALSE,3),TRUE),100))]
  LL_h <- sum(-xx_h+kk*log(xx_h)-log(factorial(kk)))

  # if(h!=0)
 (LL_h - f_res)/h
  # else
  #  y_h
}

yyy <- num_PD(y, 100, 3, h = 1e-6)

num_adj_S <- rep(0,400)
for(i in 1:400){
  num_adj_S[i] <-num_PD(y, i, 2, h = 1e-6)
}

plot(num_adj_S)
lines(adj_y[2,1,400:2], col="red")

