deriv(N) <- r * N * (N / K - 1)
deriv(adj_N) <- adj_N * r * (1- 2 * N/K)
deriv(adj_K) <- adj_N * r * (N/K)^2
deriv(adj_r) <- adj_N * N * (1 - N/K)
deriv(t_model) <- -1

initial(N) <- N_end
initial(t_model) <- t_end
initial(adj_N) <- adj_N_end
initial(adj_K) <- adj_K_end
initial(adj_r) <- adj_r_end

N_end <- user(1)
K <- user(100)
r <- user()
t_end <- user()
adj_N_end <- user()
adj_K_end <- user()
adj_r_end <- user()

