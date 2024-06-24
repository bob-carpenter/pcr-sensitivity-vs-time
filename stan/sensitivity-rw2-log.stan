data {
  int<lower = 0> N;
  array[N] int<lower = 0, upper = 1> y;
  array[N] int<lower = 0> t;
}
transformed data {
  int<lower = 2> T = max(t) + 1;
}
parameters {
  positive_ordered[T] log_theta_rev;
  real<lower = 0> sigma;
}
transformed parameters {
  vector<upper = 0>[T] log_theta = -log_theta_rev;
  vector<lower = 0, upper = 1>[T] theta = exp(log_theta);
}
model {
  sigma ~ normal(0, 0.5);
  log_theta[1] ~ normal(0, 0.05);
  log_theta[2] ~ normal(log_theta[1], 4 * sigma);
  log_theta[3:T]
    ~ normal(2 * log_theta[2:T - 1] - log_theta[1:T - 2], sigma);

  y ~ bernoulli(theta[t]);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_lpmf(y[n] | theta[t[n]]);
}
