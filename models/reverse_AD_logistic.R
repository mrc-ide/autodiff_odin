deriv(N) <- r * N * (N / K - 1)
deriv(t_model) <- -1
initial(N) <- N_end
initial(t_model) <- t_end

N_end <- user(1)
t_end <- user()
K <- user(100)
r <- user()
