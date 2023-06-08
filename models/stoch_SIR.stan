data {
  int<lower=0> T; // number of observations
  int Y[T];
  int freq;
  int<lower=0> n_steps; //number of model steps
  int n_SI[n_steps];
  int n_IR[n_steps];
}
parameters {
  real beta;
  real gamma;
}
transformed parameters {

}
model {
  int S[n_steps];
  int I[n_steps];
  int R;
  int cases_cumul;
  int cases_inc;
  real dt;
  real Y_model[T];
  real p_inf;
  real p_SI[n_steps];
  int S0;
  int I0;
  int N;

  dt = 1.0 / freq;
  S0 = 1000;

  // Initial
  S[1] = S0;
  R = 0;
  I[1] = I0;
  cases_cumul = 0;
  cases_inc = 0;

  // SIR
  for (t in 1:T) {
    for(i in 1:freq){
      int model_step = freq*(t-1)+i;
      N = S[model_step] + I[model_step] + R;
      p_inf = beta * I[model_step] / N * dt;
      p_SI[model_step] = 1 - exp(-(p_inf));

      S[model_step+1] = S[model_step] - n_SI[model_step];
      I[model_step+1] = I[model_step] + n_SI[model_step] - n_IR[model_step];
      R = R + n_IR[model_step];
      cases_cumul = cases_cumul + n_SI[model_step];
      cases_inc = cases_inc + n_SI[model_step];
    }
    Y_model[t] = cases_inc;
    cases_inc = 0;
  }

  for (t in 1:T){
    Y[t] ~ poisson(Y_model[t]);
  }
  for (t in 1:n_steps){
    n_SI[t] ~ binomial(S[t], p_SI[t]);
    n_IR[t] ~ binomial(I[t], 1 - exp(-gamma * dt));
  }

}
