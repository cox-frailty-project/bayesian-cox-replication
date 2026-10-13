data {
  int<lower=1> N, M, P;               // number of obs, participants, covariates
  vector<lower=0>[N] t;               // interval times
  array[N] int<lower=0, upper=1> C;   // Binary censored status
  array[N] int<lower=1, upper=M> id;  // patient id
  matrix[N, P] X;                     // Data matrix
}

parameters {
  // First the parameters of the baseline Weibull distribution
  real<lower=0> shape;       // Stays throughout
  real<lower=0> log_scale_0; // Updates to match individuals later
  
  // For the conclusions of interest
  vector[P] beta;       // Coefficients
  vector[M] u_unscaled; // log frailty (to be multiplied by sigma)
                        // (an attempt to reduce divergent transitions). 
  real<lower=0> sigma;  // lognormal sd, scaling factor
}

transformed parameters{
  // Scaling the log frailty
  vector[M] u = u_unscaled * sigma;

  // The Weibull scale gets updated when looking at individuals
  vector<lower=0>[N] scale_1; 
  
  // It is updated per individual as follows
  for (n in 1:N){
    scale_1[n] = exp(log_scale_0 - (X[n] * beta + u[id[n]]) / shape);
  }
}

model {
  // priors
  shape ~ normal(0,5);        // baseline Weibull shape parameter
  log_scale_0 ~ normal(0,2);  // baseline Weibull scale parameter
  beta ~ normal(0,2);         // beta coefficients
  sigma ~ normal(0,1);        // sd of lognormal distribution
  u_unscaled ~ normal(0, 1);  // sigma * u_unscaled ~ N(0, sigma^2)
  
  // observation model
  // looping over each observation
  for (n in 1:N) {
    // if it's observed, it goes to the pdf
    if (C[n] == 1) {
    target += weibull_lpdf(t[n] | shape, scale_1[n]);
    } else {
    // If it's censored it goes to the survival function
    target += weibull_lccdf(t[n] | shape, scale_1[n]);
    }
  }
}

generated quantities {
  vector[M] Z = exp(u);         // Frailties
  real sigma2 = square(sigma);  // log-frailty variance
}