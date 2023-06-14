generator <- odin::odin("models/logistic_growth_normal_obs.R")

mod <- generator$new(r=1)

tt <- seq(0, 30, length.out = 101)

tt <- seq(0, 30, length.out = 101)
y <- mod$run(tt)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")

y2 <- mod$run(tt, 50)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")
lines(y2, col = "red")
