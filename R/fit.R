library("cmdstanr")
library("ggplot2")
library("loo")
library("posterior")
library("readxl")

options(digits = 2)

printf <- function(msg, ...) { cat(sprintf(msg, ...)); cat("\n") }

df_sens <- read_excel("../data/PCR pos.xlsx")
days <- df_sens$d_sym_PCR_all
pos <- df_sens$Res_PCR
includes <- !is.na(days) & !(days < 0)
days <- days[includes]
pos <- pos[includes]
N <- length(pos)

printf("*** %d observations", N)


pos_boolean <- (pos == 1)
pos_days <- days[pos_boolean]
neg_days <- days[!pos_boolean]

df_raw <-
  data.frame(t = c(days, pos_days, neg_days),
             type = c(rep("all", length(days)),
                      rep("positive", length(pos_days)),
                      rep("negative", length(neg_days))))


plot <-
  ggplot(df_raw, aes(x = t)) +
  geom_bar(stat = "count", color = "white", size = 0.2) +
  facet_wrap(. ~ type) +
  xlab("days since symptom onset") +
  ggtitle("Counts of all tests, positive tests, and negative tests")
plot
ggsave("data-counts.pdf", plot, width = 6, height = 3)

printf("*** tests on first day = %d", sum(days == 0))
printf("*** tests on second day = %d", sum(days == 1))


printf("*** global sensitivity = %4.3f", sum(pos) / N)


mle <- function(n, N) n / N
se <- function(n, N) {
  p <- mle(n, N)
  sqrt(p * (1 - p) / N)
}

T <- max(days)
mles <- rep(NA, T)
ses <- rep(NA, T)
for (t in 0:T) {
  m <- sum(pos_days == t)
  M <- sum(days == t)
  if (M > 0) {
    mles[t + 1] <- mle(m, M)
    ses[t + 1] <- se(m, M)
  }
}

sens_breaks <- seq(0, 1, by = 0.1)
days_breaks <- seq(0, 60, by = 5)

df_mle_se <-
  data.frame(t = 0:T, mle = mles, se = ses)
plot_mle_se <-
  ggplot(df_mle_se, aes(x = t, y = mle)) +
  geom_errorbar(aes(ymin = mle - se, ymax = mle + se), color="darkgrey") +
  geom_point() +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("proportion of positive tests") +
  ggtitle("Sensitivity MLE +/- 1 std err")
plot_mle_se

ggsave("mle.pdf", plot_mle_se, width = 6, height = 3)



logit_mle <- glm(pos ~ days, family = binomial())
print(logit_mle$coefficients, digits = 2)
alpha_hat <- logit_mle$coefficients[["(Intercept)"]]
beta_hat <- logit_mle$coefficients[["days"]]


D <- max(days)
ilogit <- function(u) 1 / (1 + exp(-u))
df_log_reg <-
  data.frame(days = 0:D, sensitivity = ilogit(alpha_hat + beta_hat * (0:D)))

plot_log_reg <-
  ggplot(df_log_reg, aes(x = days, y = sensitivity)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("estimated sensitivity") +
  ggtitle("Logistic binomial regression estimate of sensitivity")
plot_log_reg


exp_mle <-
  glm(pos ~ days,
      family = binomial(link = "log"),
      start = c(-0.1, -0.1))
print(exp_mle$coefficients, digits = 2)
alpha_exp_hat <- exp_mle$coefficients[["(Intercept)"]]
beta_exp_hat <- exp_mle$coefficients[["days"]]



D <- max(days)
df_log_reg_exp <-
  data.frame(days = 0:D,
             sensitivity = exp(alpha_exp_hat + beta_exp_hat * (0:D)))

plot_log_reg_exp <-
  ggplot(df_log_reg_exp, aes(x = days, y = sensitivity)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("estimated sensitivity") +
  ggtitle("Exponential binomial regression estimate of sensitivity")
plot_log_reg_exp



printf("*** 99pct interval for ilogit(theta[d]) = (%6.4f, %6.4f)", ilogit(qnorm(0.005, 0, 3)), ilogit(qnorm(1 - 0.005, 0, 3)))


data_list <- list(N = N, y = pos, t = days + 1)
compile_fit <- function(path) {
  model <- cmdstan_model(path)
  model$sample(data = data_list, seed = 123456,
               adapt_delta = 0.95, init = 0.5, chains=4,
	       iter_warmup=10000, iter_sampling=20000, thin = 10,
               parallel_chains = 4, refresh = 5000)
}

plot_fit <- function(fit, title) {
  sens_mean <- fit$summary("theta")$mean
  sens_q5 <- fit$summary("theta")$q5
  sens_q95 <- fit$summary("theta")$q95
  days_len <- length(sens_mean)
  ggplot(data.frame(t = 0:(days_len - 1),
                    sens_mean = sens_mean,
                    sens_q5 = sens_q5,
                    sens_q95 = sens_q95),
         aes(x = t, y = sens_mean)) +
  geom_errorbar(aes(ymin = sens_q5, ymax = sens_q95), colour="darkgrey") +
  geom_point() +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("estimated sensitivity") +
  ggtitle(title)
}



sens_mix_fit <- compile_fit('../models/sensitivity/sensitivity-binomial-mix.stan')
plot_fit(sens_mix_fit,
         "Binomial mixture (logit) mean and posterior 90% interval")

sens_mix_log_fit <- compile_fit('../models/sensitivity/sensitivity-binomial-log-mix.stan')
plot_fit(sens_mix_log_fit,
         "Binomial mixture (log) mean and posterior 90% interval")

sens_shr_fit <- compile_fit('../models/sensitivity/sensitivity-shrink.stan')
plot_fit(sens_shr_fit,
         "Shrinkage model mean and posterior 90% interval")

sens_shr_mono_fit <- compile_fit("../models/sensitivity/sensitivity-shrink-mono.stan")
plot_fit(sens_shr_mono_fit,
  "Monotonic shrinkage model mean and posterior 90% interval")

sens_shr_mono_log_fit <- compile_fit("../models/sensitivity/sensitivity-shrink-mono-log.stan")
plot_fit(sens_shr_mono_log_fit,
  "Monotonic shrinkage model, log scale, mean and posterior 90% interval")

sens_rw1_fit <- compile_fit('../models/sensitivity/sensitivity-rw1.stan')
plot_fit(sens_rw1_fit,
  "Monotonic 1st-order random walk mean and posterior 90% interval")

sens_rw1_log_fit <- compile_fit('../models/sensitivity/sensitivity-rw1-log.stan')
plot_fit(sens_rw1_log_fit,
  "Monotonic 1st-order log random walk mean and posterior 90% interval")

sens_rw2_log_fit <- compile_fit('../models/sensitivity/sensitivity-rw2-log.stan')
plot_fit(sens_rw2_log_fit,
         "Monotonic 2nd-order log random walk mean and posterior 90% interval")

printf("*** posterior mean for sigma = %4.2f", mean(sens_rw1_fit$draws('sigma')))
printf("*** 90 pct central for sigma = (%4.2f, %4.2f)", quantile(sens_rw1_fit$draws('sigma'), 0.05), quantile(sens_rw1_fit$draws('sigma'), 0.95))


sens_rw2_fit <- compile_fit('../models/sensitivity/sensitivity-rw2.stan')
plot_fit(sens_rw2_fit,
  "Monotonic 2nd-order random walk mean and posterior 90% interval")
printf("*** 2nd-order RW post mean for sigma = %4.2f",
       mean(sens_rw2_fit$draws('sigma')))
printf("*** 90 pct central interval = (%4.2f, %4.2f)",
       quantile(sens_rw2_fit$draws('sigma'), 0.05),
       quantile(sens_rw2_fit$draws('sigma'), 0.95))


sens_binom_log_fit <- compile_fit('../models/sensitivity/sensitivity-binomial-log.stan')
plot_fit(sens_binom_log_fit,
  "Binomial(log) GLM mean and posterior 90% interval")
printf("*** post mean alpha binom log = %4.3f", mean(sens_binom_log_fit$draws("alpha")))
printf("*** 90 pct central = (%4.3f, %4.3f)",
       quantile(sens_binom_log_fit$draws("alpha"), 0.05),
       quantile(sens_binom_log_fit$draws("alpha"), 0.95))
printf("*** post mean beta = %4.3f", mean(sens_binom_log_fit$draws("beta")))
printf("*** central 90 pct = (%4.3f, %4.3f)",
       quantile(sens_binom_log_fit$draws("beta"), 0.05),
       quantile(sens_binom_log_fit$draws("beta"), 0.95))



sens_binom_logit_fit <- compile_fit('../models/sensitivity/sensitivity-binomial-logit.stan')
plot_fit(sens_binom_logit_fit,
  "Binomial(logit) GLM mean and posterior 90% interval")
printf("*** post mean for alpha binom logit = %4.3f",
       mean(sens_binom_logit_fit$draws("alpha")))
printf("*** central 90 pct = (%4.3f, %4.3f)",
       quantile(sens_binom_logit_fit$draws("alpha"), 0.05),
       quantile(sens_binom_logit_fit$draws("alpha"), 0.95))
printf("*** beta for binom logit = %4.3f",
       mean(sens_binom_logit_fit$draws("beta")))
printf("*** central 90 pct = (%4.3f, %4.3f)",
       quantile(sens_binom_logit_fit$draws("beta"), 0.05),
       quantile(sens_binom_logit_fit$draws("beta"), 0.95))


fit_df <- function(fit, model, link) {
  sens_mean <- fit$summary("theta")$mean
  sens_q5 <- fit$summary("theta")$q5
  sens_q95 <- fit$summary("theta")$q95
  N <- length(sens_mean)
  days <- 0:(N - 1)
  data.frame(sens_mean = sens_mean,
             sens_q5 = sens_q5,
             sens_q95 = sens_q95,
	     t = days,
	     model = rep(model, N),
	     link = rep(link, N))
}
df_all <-
  rbind(
    fit_df(sens_rw1_fit, "2. RW(1)", "logit"),
    fit_df(sens_rw2_fit, "3. RW(2)", "logit"),
    fit_df(sens_binom_logit_fit, "4. binomial", "logit"),
    fit_df(sens_shr_mono_fit, "1. heterogeneous", "logit"),
    fit_df(sens_mix_fit, "5. binomial mixture", "logit"),
    fit_df(sens_rw1_log_fit, "2. RW(1)", "log"),
    fit_df(sens_rw2_log_fit, "3. RW(2)", "log"),
    fit_df(sens_binom_log_fit, "4. binomial", "log"),
    fit_df(sens_shr_mono_log_fit, "1. heterogeneous", "log"),
    fit_df(sens_mix_log_fit, "5. binomial mixture", "log"))

comparison_plot <-
  ggplot(df_all, aes(x = t, y = sens_mean)) +
  facet_grid(link ~ model) +
  geom_line(linewidth = 0.25, color = "blue") +
  geom_errorbar(aes(ymin = sens_q5, ymax = sens_q95),
               colour="darkgrey", size = 0.25) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("estimated sensitivity")
ggsave('model-link-comparison.pdf', comparison_plot,
       width = 12,  height = 6)

df_free <-
  rbind(
    fit_df(sens_shr_fit, "shrink", "logit"),
    fit_df(sens_shr_mono_fit, "shrink (monotonic)", "logit"))
free_comparison_plot <-
  ggplot(df_free, aes(x = t, y = sens_mean)) +
  facet_grid(link ~ model) +
  geom_line(linewidth = 0.25, color = "blue") +
  geom_errorbar(aes(ymin = sens_q5, ymax = sens_q95),
               colour="darkgrey", size = 0.25) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = days_breaks) +
  scale_y_continuous(lim = c(0, 1), breaks = sens_breaks) +
  xlab("days since symptom onset") +
  ylab("estimated sensitivity")
ggsave('free-comparison.pdf', free_comparison_plot,
       width = 8,  height = 3)



# thinned to rough independence, so can ignore warning message about r_eff
draws_loo <- function(fit) {
  ll <- fit$draws("log_lik")
  loo(ll, cores = 4)
}

loo_draws <- list("binomial mixture (logit)" = draws_loo(sens_mix_fit),
       "binomial mixture (log)" = draws_loo(sens_mix_log_fit),
       "rw1 (logit)" = draws_loo(sens_rw1_fit),
       "rw1 (log)" = draws_loo(sens_rw1_log_fit),
       "rw2 (logit)" = draws_loo(sens_rw2_fit),
       "rw2 (log)" = draws_loo(sens_rw2_log_fit),
       "binomial (logit)" = draws_loo(sens_binom_logit_fit),
       "binomial (log)" = draws_loo(sens_binom_log_fit),
       "heterogeneous (logit)" = draws_loo(sens_shr_mono_fit),
       "heterogeneous (log)" = draws_loo(sens_shr_mono_log_fit))
loo_comp <- loo_compare(loo_draws)
print(loo_comp)
