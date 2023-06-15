deriv(N) <- r * N * (1 - N / K)
initial(N) <- N0

N0 <- user(1)
K <- user(100)
r <- user()

#sd_noise <- user()

#observed <- data()
#compare(observed) ~ normal(N, sd_noise)
