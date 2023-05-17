options(width = 200)
sir <- odin.dust::odin_dust("marc/sir.R")
sir_adjoint <- odin.dust::odin_dust("marc/sir_adjoint.R", debug_enable = TRUE)

read_data <- function(n, freq) {
  incidence <- read.csv("marc/incidence.csv")
  incidence <- incidence[incidence$day <= n, ]
  data <- mcstate::particle_filter_data(
    incidence, time = "day", rate = freq, initial_time = 0)
  data_input <- NULL
  for (i in seq_len(nrow(data))) {
    data_input <- c(data_input, rep(0, freq - 1L), data$cases[i])
  }
  data_input
}


f <- function(n = 100, freq = 4) {
  data <- read_data(n, freq)
  pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)
  steps <- 0:length(data)
  y <- sir$new(pars, 0, 1)$simulate(steps)[, 1, ]
  pars_adjoint <- c(pars,
                    list(main_states = y,
                         data_input = data,
                         total_steps = length(data)))
  ret <- sir_adjoint$new(pars_adjoint, 0, 1)$simulate(steps)
  ret[7:9, 1, dim(ret)[3]]
}

f()
