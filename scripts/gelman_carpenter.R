# --------------------------------------------------------------------------------
# 
#   from: https://github.com/bob-carpenter/diagnostic-testing/tree/master/src/specificity-santa-clara
# 
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(rstan)
library(nimble)
library(ggplot2)
rstan_options(auto_write = FALSE)

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
  model {
    // priors
    p ~ beta(1,1);
    spec ~ beta(1,1);
    sens ~ beta(1,1);
    // likelihood model
    real p_sample = p * sens + (1 - p) * (1 - spec);
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
  verbose = TRUE
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
  constants = list(n_sample = 3330, n_spec=401, n_sens=122)
)

sc_model_mcmc_n <- nimble::buildMCMC(conf = sc_model_n,thin = 10)

sc_model_cpp_mcmc_n <- nimble::compileNimble(sc_model_n, sc_model_mcmc_n)

sc_model_mcmc_n$run(niter = 2e4,nburnin = 500)
