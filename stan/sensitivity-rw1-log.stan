data {
  int<lower = 0> N;
  array[N] int<lower = 0, upper = 1> y;
  array[N] int<lower = 0> t;
}
transformed data {
  int<lower = 0> T = max(t) + 1;
}
parameters {
  positive_ordered[T] log_theta_rev;
  real<lower = 0> sigma;
}
transformed parameters {
  vector[T] log_theta = -log_theta_rev;
  vector[T] theta = exp(log_theta);
}
model {
  sigma ~ normal(0, 1);
  log_theta[1] ~ normal(0, 0.05);  // (0.9, 1.0) 95% interval
  log_theta[2:T] ~ normal(log_theta[1:T - 1], sigma);
  y ~ bernoulli(theta[t]);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_lpmf(y[n] | theta[t[n]]);
}
