data {
  int<lower = 0> N;
  array[N] int<lower = 0, upper = 1> y;
  array[N] int<lower = 0> t;
}
transformed data {
  int<lower = 0> T = max(t) + 1;
}
parameters {
  vector[T] logit_theta;
}
model {
  logit_theta ~ normal(0, 4);
  y ~ bernoulli_logit(logit_theta[t]);
}
generated quantities {
  vector<lower = 0, upper = 1>[T] theta
      = inv_logit(logit_theta);
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | logit_theta[t[n]]);
}
