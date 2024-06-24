data {
  int<lower = 0> N;
  array[N] int<lower = 0, upper = 1> y;
  array[N] int<lower = 0> t;
}
transformed data {
  int<lower = 2> T = max(t) + 1;
}
parameters {
  vector[2] alpha;
  ordered[2] beta;
  real<lower=0, upper=1> phi;
}
transformed parameters {
  vector<lower = 0, upper = 1>[T] theta;
  for (d in 1:T) {
    theta[d] = phi * inv_logit(alpha[1] + beta[1] * (d - 1))
      + (1 - phi) * inv_logit(alpha[2] + beta[2] * (d - 1));
  }
}
model {
  alpha ~ normal(0, 0.5);
  beta ~ normal(0, 0.5);
  phi ~ beta(2, 2);
  y ~ bernoulli(theta[t]);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(y[n] | theta[t[n]]);
  }
}
