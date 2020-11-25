# --------------------------------------------------------------------------------
# 
#   from: https://github.com/bob-carpenter/diagnostic-testing/tree/master/src/specificity-santa-clara
# 
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()
library(rstan)
library(nimble)
library(ggplot2)
rstan_options(auto_write = FALSE)


# --------------------------------------------------------------------------------
#   simplest model, no poststratification
# --------------------------------------------------------------------------------

santa_clara <- '
  data {
    int<lower = 0> y_sample;
    int<lower = 0> n_sample;
    int<lower = 0> y_spec;
    int<lower = 0> n_spec;
    int<lower = 0> y_sens;
    int<lower = 0> n_sens;
  }
  parameters {
    real<lower=0, upper = 1> p;
    real<lower=0, upper = 1> spec;
    real<lower=0, upper = 1> sens;
  }
  transformed parameters {
    real<lower=0, upper=1> p_sample = p * sens + (1 - p) * (1 - spec);
  }
  model {
    p ~ beta(1,1);
    spec ~ beta(1,1);
    sens ~ beta(1,1);
    y_sample ~ binomial(n_sample, p_sample);
    y_spec ~ binomial(n_spec, spec);
    y_sens ~ binomial(n_sens, sens);
  }
'

sc_data <- list(y_sample=50, n_sample=3330, y_spec=369+30, n_spec=371+30, y_sens=25+78, n_sens=37+85)

sc_model <- rstan::stan(
  model_code = santa_clara,
  data = sc_data,
  chains = 4,
  warmup = 1e4,
  iter = 2e4,
  cores = 4,
  verbose = F
)

stan_trace(sc_model)
stan_dens(sc_model)
stan_par(sc_model,"sample")

draws_1 <- extract(sc_model)
x <- as.vector(draws_1$spec)
y <- as.vector(draws_1$p)
hist(y)

# par(mar=c(3,3,0,1), mgp=c(2, .7, 0), tck=-.02)
# plot(x, y, xlim=c(min(x), 1), ylim=c(0, max(y)), xaxs="i", yaxs="i", xlab=expression(paste("Specificity, ", gamma)), ylab=expression(paste("Prevalence, ", pi)), bty="l", pch=20, cex=.3)
# dev.off()



santa_clara_n <- nimble::nimbleCode({
  p ~ dbeta(shape1 = 1, shape2 = 1)
  spec ~ dbeta(shape1 = 1, shape2 = 1)
  sens ~ dbeta(shape1 = 1, shape2 = 1)
  p_sample <- p * sens + (1 - p) * (1 - spec)
  y_sample ~ dbinom(prob = p_sample, size = n_sample)
  y_spec ~ dbinom(prob = spec, size = n_spec)
  y_sens ~ dbinom(prob = sens, size = n_sens)
})


sc_model_n <- nimble::nimbleModel(
  code = santa_clara_n,
  data = list(y_sample=50, y_spec=369+30, y_sens=25+78),
  constants = list(n_sample = 3330, n_spec=401, n_sens=122),
  inits = list(p = 0.1, spec = 0.9, sens= 0.9)
)

sc_model_mcmc_n <- nimble::buildMCMC(conf = sc_model_n,thin = 10)

sc_model_cpp_n <- nimble::compileNimble(sc_model_n)

sc_model_cpp_mcmc_n <- nimble::compileNimble(sc_model_mcmc_n, project = sc_model_n)

samples <- nimble::runMCMC(mcmc = sc_model_cpp_mcmc_n,niter = 1e5,nburnin = 1e4,thin = 10,progressBar = TRUE)

hist(samples[,"p"])


# --------------------------------------------------------------------------------
#   better model: hierarchical
# --------------------------------------------------------------------------------


sc_hier <- '
  data {
    int<lower = 0> y_sample;
    int<lower = 0> n_sample;
    int<lower = 0> J_spec;
    int<lower = 0> y_spec[J_spec];
    int<lower = 0> n_spec[J_spec];
    int<lower = 0> J_sens;
    int<lower = 0> y_sens[J_sens];
    int<lower = 0> n_sens[J_sens];
    real<lower = 0> logit_spec_prior_scale;
    real<lower = 0> logit_sens_prior_scale;
  }
  parameters {
    real<lower = 0, upper = 1> p;
    real mu_logit_spec;
    real mu_logit_sens;
    real<lower = 0> sigma_logit_spec;
    real<lower = 0> sigma_logit_sens;
    vector<offset = mu_logit_spec, multiplier = sigma_logit_spec>[J_spec] logit_spec; // logit(gamma_j)
    vector<offset = mu_logit_sens, multiplier = sigma_logit_sens>[J_sens] logit_sens; // logit(delta_j)
  }
  transformed parameters {
    vector[J_spec] spec = inv_logit(logit_spec);
    vector[J_sens] sens = inv_logit(logit_sens);
  }
  model {
    real p_sample = p * sens[1] + (1 - p) * (1 - spec[1]); // SC frequency of positive tests (prev p corrected for imperfect tests)
    y_sample ~ binomial(n_sample, p_sample); // SC data likelihood
    y_spec ~ binomial(n_spec, spec); // likelihood for 13 spec studies
    y_sens ~ binomial(n_sens, sens); // likelihood for 3 sens studies
    logit_spec ~ normal(mu_logit_spec, sigma_logit_spec); // spec prior is normal on log-odds scale
    logit_sens ~ normal(mu_logit_sens, sigma_logit_sens); // sens prior is normal on log-odds scale
    sigma_logit_spec ~ normal(0, logit_spec_prior_scale); // hyperprior
    sigma_logit_sens ~ normal(0, logit_sens_prior_scale); // hyperprior
    mu_logit_spec ~ normal(4, 2);  // weak hyperprior on mean of distribution of spec
    mu_logit_sens ~ normal(4, 2);  // weak hyperprior on mean of distribution of sens
  }
'

sc_hier_data <- list(
  # SC data
  y_sample=50,
  n_sample=3330,
  # data from 13 specificity studies in Bendavid 2020b
  J_spec=14,
  y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
  n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
  # data from 3 sensitivity studies in Bendavid 2020b
  J_sens=4,
  y_sens=c(0, 78, 27, 25),
  n_sens=c(0, 85, 37, 35)
)

# Fit with weak priors
sc_hier_data$logit_spec_prior_scale <- 1
sc_hier_data$logit_sens_prior_scale <- 1

sc_hier_model <- rstan::stan(
  model_code = sc_hier,
  data = sc_hier_data,
  chains = 4,
  warmup = 1e4,
  iter = 2e4,
  cores = 4,
  control=list(adapt_delta=0.95),
  verbose = TRUE
)

stan_dens(sc_hier_model)
