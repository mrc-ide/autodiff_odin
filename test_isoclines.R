#Mean of the two distributions
mu <- c(prior_par$meanlog_beta,prior_par$meanlog_gamma)
#Build VCV matrix assuming independance between the two prior distributions
l <- matrix(c(prior_par$sdlog_beta^2,0,0,prior_par$sdlog_gamma^2), nrow=2)
#Draw sample from prior
draw_prior <- mvtnorm::rmvnorm(n = 1000, mean = mu, sigma = l)

posterior_value_samples <- calculate_posterior_map(draw_prior,
                                                   mcmc_pars,
                                                   filter)

n_pal <- 20
rbPal <- colorRampPalette(c('yellow','blue'))
plot(draw_prior, pch = 19, col=grey(.8),
     xlab="log(beta)", ylab="log(gamma)", # xlim = c(-5,3), ylim = c(-5,3),
     main ="Variational inference")
legend("bottom",title="Ventile of LogPosterior",
       legend=c(1:n_pal),col =rbPal(n_pal),pch=20, horiz = TRUE, bty = "n", cex = .6)
#plot points with the colours
finite_prior <- posterior_value_samples!=-Inf
colZ <- rbPal(n_pal)[as.numeric(
  cut(log(abs(posterior_value_samples[finite_prior])),breaks = n_pal))]
points(draw_prior[finite_prior,],pch = 20,col = colZ)


#Plot 95% CI of the prior distribution
lines(ellipse::ellipse( l , centre = mu) , col='red', lwd=3, lty = 3)

#Plot the mcmc samples
points(log(as.numeric(mcmc_run$pars[,"beta"])),
       log(as.numeric(mcmc_run$pars[,"gamma"])),
       xlim=c(-7,1), ylim = c(-5,1), col="red", pch=19)

iter_grad <- function(x, y, n_step, l){
  for(i in seq(n_step)){
    grad_theta <- gradient_LP(c(x,y), mcmc_pars, filter)
    g_x <- grad_theta$grad_LP[1]
    g_y <- grad_theta$grad_LP[2]
    x_h <- x + l*g_y/sqrt(g_x^2+g_y^2)
    y_h <- y - l*g_x/sqrt(g_x^2+g_y^2)
    lines(c(x, x_h),c(y,y_h))
    x<-x_h
    y<-y_h
  }
}

#Draw one random point from the prior
n_iso <- 10
sample_n <- sample(nrow(draw_prior),n_iso)
#sample_n <- which.max(posterior_value_samples)

points(draw_prior[sample_n,1], draw_prior[sample_n,2], col="orange", pch=19)

l <- .05
n_step <- 500

for(i in sample_n)
{
  x <- draw_prior[i,1]
  y <- draw_prior[i,2]
  iter_grad(x,y,n_step,l)
  iter_grad(x,y,n_step,-l)
}




