p_IR <- 1 - exp(-(gamma))
S0 <- 1000
freq <- user(4)
dt <- 1.0 / freq

N <- S + I + R
p_inf <- beta * I / N
p_SI <- 1 - exp(-(p_inf))
n_SI <- S * p_SI * dt
n_IR <- I * p_IR * dt

update(time) <- (step + 1) * dt
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(cases_cumul) <- cases_cumul + n_SI
update(cases_inc) <- if (step %% freq == 0) n_SI else cases_inc + n_SI

initial(time) <- 0
initial(S) <- S0
initial(R) <- 0
initial(I) <- I0
initial(cases_cumul) <- 0
initial(cases_inc) <- 0

beta <- user(0.2)
gamma <- user(0.1)
I0 <- user(10)
