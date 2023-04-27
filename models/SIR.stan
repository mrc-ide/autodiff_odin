data {
  int<lower=0> T; // number of observations
  int Y[T];
  int freq;
}
parameters {
  real beta;
  real gamma;
  real I0;
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

  dt = 1.0 / freq;
  p_IR = 1 - exp(-(gamma * dt));
  S0 = 1000;

  // Initial
  S = S0;
  R = 0;
  I = I0;
  cases_cumul = 0;
  cases_inc = 0;

  // SIR
  for (t in 1:T) {
    for(i in 1:freq){
      N = S + I + R;
      p_inf = beta * I / N * dt;
      p_SI = 1 - exp(-(p_inf));
      n_SI = S * p_SI;
      n_IR = I * p_IR;

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

}
