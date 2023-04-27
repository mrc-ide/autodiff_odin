incidence <- read.csv("data/incidence.csv")
data <- mcstate::particle_filter_data(
  incidence, time = "day", rate = 4, initial_time = 0)

sir_stoch <- odin.dust::odin_dust("models/sir_stochastic.R")

beta <- 0.25
gamma <- 0.1
I0 <- 1
pars <- list(beta = beta, gamma = gamma, I0 = I0)
sir_stoch_mod <- sir_stoch$new(pars, 0, 100)

y <- sir_stoch_mod$simulate(c(0, data$time_end))
i <- sir_stoch_mod$info()$index[["time"]]
j <- sir_stoch_mod$info()$index[["cases_inc"]]
matplot(y[i, 1, ], t(y[j, , ]), type = "l", col = "#00000055", lty = 1, las = 1,
        xlab = "Day", ylab = "Cases")
