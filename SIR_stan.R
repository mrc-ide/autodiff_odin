model_SIR <-
'data {
  int<lower=0> T; // number of observations
  int Y[T];
  int freq;
}
parameters {
  real<lower=0,upper=4> beta;
  real<lower=0,upper=5> gamma;
  real<lower=0,upper=100> I0;
}
transformed parameters {
  real S;
  real I;
  real R;
  real cases_cumul;
  real cases_inc;
  real N;
  real p_IR;
  real S0;
  real dt;
  real Y_model[T];
  real p_inf;
  real p_SI;
  real n_SI;
  real n_IR;

  p_IR = 1 - exp(-(gamma));
  S0 = 1000;
  dt = 1.0 / freq;

  // Initial
  S = S0;
  R = 0;
  I = 0;
  cases_cumul = 0;
  cases_inc = 0;

  // SIR
  for (t in 1:T) {
    for(i in 1:freq){
      N = S + I + R;
      p_inf = beta * I / N;
      p_SI = 1 - exp(-(p_inf));
      n_SI = S * p_SI * dt;
      n_IR = I * p_IR * dt;

      S = S - n_SI;
      I = I + n_SI - n_IR;
      R = R + n_IR;
      cases_cumul = cases_cumul + n_SI;
      cases_inc = cases_inc + n_SI;
    }
    Y_model[t] = cases_inc;
    cases_inc = 0;
  }
}
model {

  for (t in 1:T){
    Y[t] ~ poisson(Y_model[t]);
  }
}'

m <- stan_model(model_code = model_SIR)

opt <- optimizing(m, data = data_lst)
init_lst <- purrr::map(1:4, function(i) {
  list(
    beta = rnorm(3,opt$par[1:3], 0.01),
    gamma = opt$par[4],
    thetap = opt$par[5:7],
    thetan = opt$par[8],
    I0 = opt$par[9:11],
    susc = opt$par[12:14],
    asc = opt$par[15],
    niliPeriod = opt$par[16],
    niliAmplitude = opt$par[17],
    niliPeak = opt$par[18]
  )})
stan_fit <- rstan::sampling(m, data = data_lst, chains = 4, iter = 2000, thin = 1, cores = 4, control = list(adapt_delta = 0.95, max_treedepth = 15), init = init_lst)
multistrain_fit <- rstan::extract(stan_fit)
