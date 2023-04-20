model_SIR <-
'data {
  int<lower=0> N; // number of observations
  int Y[N];
}
parameters {
  real<lower=0,upper=4> beta;
  real<lower=0,upper=5> gamma;
  rel<lower=0,upper=100> I0;
}
transformed parameters {
  real S[K,T];
  real I[K,T];
  real R[K,T];
  matrix[K,W] PILI;
  real NILI[W];
  real ILI[W];
  matrix[K,T-1] pILI;
  vector<lower=0>[T-1] nILI;
  real<lower=0> newI[K,T-1];
  real<lower=0> newR[K,T-1];

  vector[K+1] pPILI[W];

  for (w in 1:W) NILI[w] = 0;

  // Initial
  for (k in 1:K) {
    I[k,1] = exp(I0[k])*susc[k]*N;
    S[k,1] = susc[k]*N - I[k,1];
    R[k,1] = (1-susc[k])*N;
    for (w in 1:W) PILI[k,w] = 0;
  }

  // SIR
  for (t in 2:T) {
    for (k in 1:K) {
      newI[k,t-1] = fmin(S[k,t-1], S[k,t-1]*beta[k]*I[k,t-1]/N);

      S[k,t] = S[k,t-1] - newI[k,t-1];
      I[k,t] = I[k,t-1] + newI[k,t-1];

      newR[k,t-1] = fmin(I[k,t], I[k,t]*gamma);
      I[k,t] = I[k,t] - newR[k,t-1];
      R[k,t] = R[k,t-1] + newR[k,t-1];


      pILI[k,t-1] = newI[k,t-1]*thetap[k];
    }
    // Uninfected multiplied by the influenza negative ILI rate
    nILI[t-1] = fmax(0,N-sum(I[,t]))*exp(fmin(thetan + niliForcing(t, niliAmplitude, niliPeriod, niliPeak), 0));
  }

  // Per week
  for (w in 1:W) {
    for (t in ((w-1)*7+1):(w*7)) {
      for (k in 1:K) {
        // Calculate pili and nili
        PILI[k,w] = PILI[k,w] + pILI[k,t];
      }
      NILI[w] = NILI[w] + nILI[t];
    }
    ILI[w] = fmax(sum(col(PILI, w)) + NILI[w],1);
    for (k in 1:K) pPILI[w,k] = PILI[k,w]/ILI[w];
    pPILI[w,K+1] = 1 - sum(pPILI[w,1:K]);
  }
}
model {
  real sheddingPeriod;
  sheddingPeriod = 1/gamma;
  sheddingPeriod ~ normal(4.8, 0.245); // carrat_time_2008


  for (k in 1:K) {
    real reff;
    reff = beta[k]*susc[k]/gamma;
    reff ~ normal(1.28, 0.133); // Assume biggerstaff is on Reff
  }

  thetap[1] ~ beta(18.21082, 30.61019); // H1N1
  thetap[2] ~ beta(36.41735, 52.98367); // H3N2
  thetap[3] ~ beta(4.54816, 50.90341); // B

  asc ~ beta(35.644, 69.314);

  for (w in 1:W) {
    YILI[w,1] ~ approxbin(YILI[w,2]*ILI[w]/N, asc);
    Y[w] ~ multinomial(pPILI[w]);
  }
}'

m <- stan_model(model_code = model_str)

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
