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
}
transformed parameters {
  vector<upper = 0>[T] log_theta = -log_theta_rev;
  vector<lower = 0, upper = 1>[T] theta = exp(log_theta);
}
model {
  log_theta ~ normal(0, 3);
  y ~ bernoulli(theta[t]);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_lpmf(y[n] | theta[t[n]]);
}
