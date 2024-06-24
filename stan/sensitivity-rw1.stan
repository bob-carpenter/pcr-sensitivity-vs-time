data {
  int<lower = 0> N;
  array[N] int<lower = 0, upper = 1> y;
  array[N] int<lower = 0> t;
}
transformed data {
  int<lower = 0> T = max(t) + 1;
}
parameters {
  ordered[T] logit_theta_rev;
  real<lower = 0> sigma;
}
transformed parameters {
  vector[T] logit_theta = -logit_theta_rev;
}
model {
  sigma ~ normal(0, 1);
  logit_theta[1] ~ normal(4, 1);
  logit_theta[2:T] ~ normal(logit_theta[1:T - 1], sigma);
  y ~ bernoulli_logit(logit_theta[t]);
}
generated quantities {
  vector<lower = 0, upper = 1>[T] theta
      = inv_logit(logit_theta);
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit_theta[t[n]]);
}
