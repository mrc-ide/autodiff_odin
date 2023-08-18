# remotes::install_github("mrc-ide/odin@mrc-4358", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/odin.dust@mrc-4359", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/dust@mrc-4307", upgrade = FALSE, force = TRUE)
# remotes::install_github("mrc-ide/mcstate", upgrade = FALSE, force = TRUE)

# Create generator for dust model
gen <- odin.dust::odin_dust("models/sir_adjoint.R")

# Create dust data
incidence <- read.csv("data/incidence.csv")
incidence <- data.frame(
  time = incidence$day * 4,
  cases_observed = incidence$cases)
d <- dust::dust_data(incidence)

# Create set of parameters
pars <- list(beta = 0.25, gamma = 0.1, I0 = 1)

# Create new deterministic model with data attached
mod <- gen$new(pars, 0, 1, deterministic = TRUE)
mod$set_data(d)

## This is the current temporary arrangement with dust and may change:
info <- mod$info()

compute_gradient <- function(mod, pars, t = 0){
  mod$update_state(pars, time = t)
  mod$run_adjoint()
}

n_samples <- 1000
fixed_I0 <- 5
prior_df <- data.frame(beta = rlnorm(n_samples, meanlog = -1, sdlog = 0.4),
                       gamma = rlnorm(n_samples, meanlog = -1, sdlog = 0.4),
                       I0 = rep(fixed_I0, n_samples)) #rlnorm(n_samples))

df <- NULL
for(i in seq(nrow(prior_df))){
  res <- compute_gradient(mod,as.list(prior_df[i,]))
  if(!is.infinite(res$log_likelihood))
    df <- rbind(df,
                data.frame(cbind(prior_df[i,],
                                 log_likelihood = res$log_likelihood,
                                 gradient_beta = res$gradient[1],
                                 gradient_gamma = res$gradient[2],
                                 gradient_I0 = res$gradient[3])))
}
row.names(df) <- NULL

LL_5q <- quantile(df$log_likelihood,0.95)
good_LL <- df$log_likelihood > LL_5q
plot(prior_df$beta, prior_df$gamma,
     pch = 19, col=grey(.9), xlab="beta", ylab="gamma",
     xlim = c(0,1.5), ylim = c(0,1.5))
points(df$beta, df$gamma, pch = 19, col=grey(.4))
points(df$beta[good_LL], df$gamma[good_LL], pch = 19, col="red")

up_target <- min(df$log_likelihood[good_LL])
low_target <- max(df$log_likelihood[!good_LL])
mid_point <- which(df$log_likelihood==low_target)

points(df[mid_point,"beta"], df[mid_point,"gamma"], pch=19, col="orange")

theta <- df[mid_point,c("beta","gamma","I0")]
dist <- function(x) sqrt(sum(x["beta"]^2+x["gamma"]^2))

l <- 0.0001
n_step <- 10000
LP_chain <- NULL
LL_chain <- NULL
for(i in seq(n_step)){
  grad_theta <- compute_gradient(mod,as.list(theta))

  theta_h <- theta + l*c(grad_theta$gradient["gamma"],
                         -grad_theta$gradient["beta"],
                         0)/dist(grad_theta$gradient)
  lines(c(theta["beta"], theta_h["beta"]), c(theta["gamma"], theta_h["gamma"]), col="green", lwd = 3)

  theta <- theta_h

  LP_chain <- rbind(LP_chain,theta)
  LL_chain <- rbind(LL_chain, grad_theta$log_likelihood)

}

#plot(seq(from=l, to=l*n_step, by=l),LL_chain, type="l")
#lines(seq(from=l, to=l*n_step, by=l),LL_chain, col="blue")











# x <- c(0.25,0.1,1)
# h <- 1e-7
# names(x) <- c("beta","gamma","I0")
# LL_chain <- NULL
#
# for(i in 1:1000){
#   res <- mod$run_adjoint()
#   LL_chain <- c(LL_chain, res$log_likelihood)
#   x <- x + h * res$gradient
#   mod$update_state(list(beta = x["beta"], gamma = x["gamma"], I0 = x["I0"]), time = 0)
# }
#
# x <- c(0.25, 0.1, 1)
# h <- 1e-8
# names(x) <- c("beta", "gamma", "I0")
# LL_chain2 <- NULL
#
# # Initialize NAG-specific variables
# momentum <- 0.8
# velocity <- rep(0, length(x))
#
# for (i in 1:1000) {
#   # Calculate gradient using the current lookahead position
#
#   lookahead_x <- x + momentum * velocity
#   mod$update_state(list(beta = lookahead_x["beta"], gamma = lookahead_x["gamma"], I0 = lookahead_x["I0"]), time = 0)
#   res <- mod$run_adjoint()
#   gradient <- res$gradient
#
#   # Update Nesterov velocity
#   velocity <- momentum * velocity + h * gradient
#
#   # Update the parameters using the Nesterov velocity
#   x <- x + velocity
#
#   # Update model state
#   mod$update_state(list(beta = x["beta"], gamma = x["gamma"], I0 = x["I0"]), time = 0)
#   res <- mod$run_adjoint()
#   LL_chain2 <- c(LL_chain2, res$log_likelihood)
# }
