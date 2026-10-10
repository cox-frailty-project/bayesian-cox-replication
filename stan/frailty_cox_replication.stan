// Bayesian Cox model with shared Gaussian frailty (McGilchrist & Aisbett, 1991)
// Requires data sorted by ascending time t.
data {
  int<lower=1> N, M, P;                  // intervals, patients, covariates
  vector[N] t;                           // gap times, sorted ascending
  array[N] int<lower=0, upper=1> d;      // 1 = infection, 0 = censored
  array[N] int<lower=1, upper=M> pid;    // patient_id for each interval
  matrix[N, P] X;                        // covariates- age,sex,disease type
  real<lower=0> tau;                     // fixed prior sd of beta
}
transformed data {
  // Breslow ties: Tied events are all compared against the same set of patients still at risk.
  array[N] int first;                    // first row of each tie block
  first[1] = 1;
  for (k in 2:N)
    first[k] = (t[k] == t[k - 1]) ? first[k - 1] : k;
}
parameters {
  vector[P] beta;                        // regression coeffs
  sum_to_zero_vector[M] u;               // log-frailties restriction [1'u = 0]
  real<lower=0> sigma;                   // log-frailty sd
}
transformed parameters {
  vector[N] eta = X * beta + u[pid];     // linear substitution
}
model {
  // Priors
  beta ~ normal(0, tau);                 // fixed-tau ridge prior
  sigma ~ normal(0, 1);                  // half-normal via lower = 0
  u ~ normal(0, sigma);                  // log-frailties ~ normal(0, sigma^2)
  target += log(sigma);                  // exponent M -> M-1, since 1'u = 0

  // Cox log-partial likelihood (Breslow ties):
  for (k in 1:N)
    if (d[k] == 1)
      target += eta[k] - log_sum_exp(eta[first[k]:N]);
}
generated quantities {
  vector[M] Z = exp(u);                  // Frailties
  real sigma2 = square(sigma);           // log-frailty sd
}

