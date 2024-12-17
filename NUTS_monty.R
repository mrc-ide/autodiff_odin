# Load model
m <- monty::monty_example("banana", sigma = 0.5)

# Check value of gradient
monty::monty_model_gradient(m, c(1,-1))

# Compute the hamiltonian i.e. kinetic + potential energies
hamiltonian <- function(model, theta, r){
  #Note that we only need the likelihood here so no need to also calculate
  #the gradient as we do here
  sum(r^2) / 2  - monty::monty_model_density(model, theta)
}

# Check Hamiltonian value
hamiltonian(m, c(1,-1), c(1,-1))

# Perform 1 leafrog integration of step epsilon
leapfrog <- function(model, current_theta, current_r, epsilon){
  # initialise to the current value of theta and r
  theta <- current_theta
  r <- current_r
  # Make a half step for momentum
  r <- r + epsilon * monty::monty_model_gradient(m, theta) / 2
  # Make a full step for theta
  theta <- theta + epsilon * r
  # Make a half step for momentum
  r <- r + epsilon * monty::monty_model_gradient(m, theta) / 2
  return(list(theta = theta, r = r))
}

# Recursive construction of sample exploration tree
build_tree <- function(model, theta, r, u, v, j, epsilon, theta_0, r_0, delta = 1000)
{
  if(j==0) {
    # base case, one leapfrog in direction v
    theta_r_prop <- leapfrog(model, theta, r, v*epsilon)
    H_prop <- hamiltonian(model, theta_r_prop$theta, theta_r_prop$r)
    H_0 <- hamiltonian(model, theta_0, r_0)
    n <- as.integer(u <= exp(-H_prop))
    #lines(c(theta["beta"],theta_r_prop$theta["beta"]),c(theta["gamma"],theta_r_prop$theta["gamma"]),col=grey(.6))
    #points(theta_r_prop$theta["beta"],theta_r_prop$theta["gamma"], col=grey(.6))
    s <- u < exp(delta - H_prop)
    return(list(theta_minus = theta_r_prop$theta,
                r_minus = theta_r_prop$r,
                theta_plus = theta_r_prop$theta,
                r_plus = theta_r_prop$r,
                theta_prop = theta_r_prop$theta,
                n_prop = n,
                s_prop = s,
                alpha = min(1, exp(H_0-H_prop)),
                n_alpha = 1)
    )
  } else { #j>0
    result_list <- build_tree(model, theta, r, u, v, j-1, epsilon, theta_0, r_0, delta)
    if(result_list$s_prop){ # continue the tree unless stop condition is reached
      if(v==-1){
        alternative_list <- build_tree(model, result_list$theta_minus, result_list$r_minus,
                                       u, v, j-1, epsilon, theta_0, r_0, delta)
        result_list$theta_minus <- alternative_list$theta_minus
        result_list$r_minus <- alternative_list$r_minus
      } else { #v==1
        alternative_list <- build_tree(model, result_list$theta_plus, result_list$r_plus,
                                       u, v, j-1, epsilon, theta_0, r_0, delta)
        result_list$theta_plus <- alternative_list$theta_plus
        result_list$r_plus <- alternative_list$r_plus
      }
      sum_n_prop <- result_list$n_prop+alternative_list$n_prop
      if(sum_n_prop > 0)
        if(runif(1)<alternative_list$n_prop/sum_n_prop)
          result_list$theta_prop <- alternative_list$theta_prop
      result_list$alpha <- result_list$alpha + alternative_list$alpha
      result_list$n_alpha <- result_list$n_alpha + alternative_list$n_alpha
      result_list$s_prop <- alternative_list$s_prop &
        ((result_list$theta_plus-result_list$theta_minus)%*%result_list$r_minus >= 0) &
        ((result_list$theta_plus-result_list$theta_minus)%*%result_list$r_plus >= 0)
      result_list$n_prop <- sum_n_prop
    }
    return(result_list)
  }
}

# NUTS step
NUTS_step <- function(model, theta, epsilon, D_max){
  theta_prop <- theta
  r0 <- rnorm(length(theta),0,1)
  u <- runif(1)*exp(-hamiltonian(model, theta, r0))
  tree_list <- list(theta_minus = theta,
                    r_minus = r0,
                    theta_plus = theta,
                    r_plus = r0
  )
  j <- 0
  n <- 1
  s <- TRUE
  while(s){
    v <- sample(c(-1,1),1)
    if(v==-1){
      tree_list <- build_tree(model, tree_list$theta_minus, tree_list$r_minus,
                              u, v, j, epsilon, theta, r0, D_max)
    } else {
      tree_list <- build_tree(model, tree_list$theta_plus, tree_list$r_plus,
                              u, v, j, epsilon, theta, r0, D_max)
    }
    if(tree_list$s_prop)
      if(runif(1)<min(1,tree_list$n_prop/n))
        theta_prop <- tree_list$theta_prop
    n <- n + tree_list$n_prop
    s <- tree_list$s_prop &
      ((tree_list$theta_plus-tree_list$theta_minus)%*%tree_list$r_minus >= 0) &
      ((tree_list$theta_plus-tree_list$theta_minus)%*%tree_list$r_plus >= 0)
    j <- j+1
  }
  list(theta_prop=theta_prop, j=j, s=s, n=n, theta_minus=tree_list$theta_minus, theta_plus=tree_list$theta_plus)
}

# Sampler loop
theta0 <- c(1, -1)
M <- 5000
D_max <- 1000
#epsilon0 <- find_epsilon1(mod, theta0, g, dg, 0.0001)
epsilon0 <- 0.018
mu <- log(10*epsilon0)/10
theta_m <- matrix(rep(theta0, M+1), ncol = length(theta0), byrow = TRUE)
colnames(theta_m) <- names(theta0)

for(i in 1:M)
{
  res <- NUTS_step(m, theta_m[i,], epsilon0, D_max)
  theta_m[i+1,] <- res$theta_prop
}

# Visualisation of samples
a <- seq(-2, 6, length.out = 1000)
b <- seq(-2, 2, length.out = 1000)
z <- outer(a, b, function(alpha, beta) {
  exp(monty::monty_model_density(m, rbind(alpha, beta)))
})

theta <- seq(0, 2 * pi, length.out = 10000)
z95 <- local({
  sigma <- 0.5
  r <- sqrt(qchisq(.95, df = 2))
  x <- r * cos(theta)
  y <- r * sin(theta)
  cbind(x^2 + y * sigma, x)
})
image(a, b, z, xlab = "alpha", ylab = "beta")
points(theta_m[-(1:50),1], theta_m[-(1:50),2],pch=19, col="#3322ff44")
lines(z95[, 1], z95[, 2])
