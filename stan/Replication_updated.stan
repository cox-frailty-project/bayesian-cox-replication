data {
  int<lower=1> M, N, P, m; // # of patients, total obs, # of covariates, # of 
  // recurrences per patient
  array[N] int<lower=0, upper=1> C; // C = 1 means observed, C = 0 means censored
  array[N] int<lower=1, upper=M> pid; // patient_id for each interval
  matrix[N, P] X; // covariates- age,sex,disease type
  matrix[N,M] r; // We introduce this r matrix for calculating alphas (r as in the paper)
  array[M,m] int a_ind; // Also for calculating alphas (what rows the occurrences happen)
  real<lower=0> sigma; // The sd, we initialise it ourselves (one option)
}
parameters {
  vector[P] beta; // regression coeffs
  sum_to_zero_vector[M] u; // log-frailties restriction [1'u = 0]
  // real<lower=0> sigma; // If we want sigma as a parameter instead        
}
transformed parameters {
  // We throw in alpha here.
  // alpha is defined at every time point tn and every individual (pid)
  // So it's an N x M matrix (or MxN depending on what order you want it).
  
  // I think it will help if we do a couple things. 1) Define r_{i,n} (done),
  // and 2) Define g(x,b) values in a matrix
  
  // We start with the g(x,beta) values for the individuals at their event times
  matrix[M,m] g;
  for (k in 1:M){
    // a_ind is defined in the data, it says where the event times are for 
    // person i in the k dataset (and correspondingly the X covariates)
    g[k,1] = exp(X[a_ind[k,1],]*beta);
    g[k,2] = exp(X[a_ind[k,2],]*beta);
  } 
  
  // Then the next big thing is getting alpha
  matrix[N,M] alpha;
  // We loop over the observations (rows)
  for (i in 1:N){
    // And fill each column at each row
    for (k in 1:M){
      // And honestly an if else has to be the easiest option
      if (r[i,k] == 1){
        alpha[i,k] = g[k,1] + g[k,2];
      } else if (r[i,k] == 2){
        alpha[i,k] = g[k,2];
      } else {
        alpha[i,k] = 0.0;
      }
    }
  } // Now we should have a matrix of alphas

  // Now we want the b's. This should also be a matrix NxM I'd say.
  
  // I think it's nice to do the denominators of that expression first
  // i.e. the total hazard of everyone left at tn. So every time point you have 
  // one denom, so we're summing over the columns at every tn.
  vector[N] b_denom;
  matrix[N,M] b_denom_parts;
  for (i in 1:N){
    for (k in 1:M){
      b_denom_parts[i,k] = exp(u[k])*alpha[i,k];
    }
  }
  for (i in 1:N){
    b_denom[i] = sum(b_denom_parts[i,]);
  } 
  
  // And now calculating them all proper
  matrix[N,M] b;
  for (i in 1:N){
    for (k in 1:M){
      b[i,k] = alpha[i,k]/b_denom[i];
    }
  } // YAY
}
model {
  // Priors
  beta ~ normal(0, 10);  // for the coefs
  // sigma ~ normal(0, 50); // A prior for sigma if we had it as a parameter
  // u ~ normal(0, sigma); // A prior for the frailties
  target += -0.5 * dot_self(u) / square(sigma); // This is the penalty part
  // of the likelihood. If we define it this way, we don't include the prior
  // for the log frailties or sigma. 

  // THE LIKELIHOOD
  for (i in 1:N)
    if (C[i] == 1)
      // Fixed form
      target += u[pid[i]] + log(b[i,pid[i]]);
}
generated quantities {
  vector[M] Z = exp(u);                  // Frailties
  //real sigma2 = square(sigma);         // log-frailty sd
  real sigma2 = dot_self(u)/(M-1); // We calculate the variance as the variance 
  // of the frailties, not as the square of our old parameter.
}
