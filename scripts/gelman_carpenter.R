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

# STAN
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

# NIMBLE
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
    p ~ beta(1,1);
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

# NIMBLE
sc_hier_n <- nimble::nimbleCode({
  
  # prior models
  
  # priors (scalar)
  p ~ dbeta(shape1 = 1, shape2 = 1) # uniform prior on prevalence
  mu_logit_spec ~ dnorm(mean = 4, sd = 2) # hyperprior on mean of distribution of spec (gamma)
  mu_logit_sens ~ dnorm(mean = 4, sd = 2) # hyperprior on mean of distribution of sens (delta)
  sigma_logit_spec ~ T(dnorm(mean = 0, sd = 1),min = 0) # hyperprior on sd of distribution of spec (gamma)
  sigma_logit_sens ~ T(dnorm(mean = 0, sd = 1),min = 0) # hyperprior on sd of distribution of sens (delta)
  # priors (for studies)
  for(i in 1:J_spec){
    logit_spec[i] ~ dnorm(mean = mu_logit_spec, sd = sigma_logit_spec) # spec prior is normal on log-odds scale
  }
  for(i in 1:J_sens){
    logit_sens[i] ~ dnorm(mean = mu_logit_sens, sd = sigma_logit_sens) # sens prior is normal on log-odds scale
  }
  
  # likelihood models
  
  # likelihood for 13 spec studies
  for(i in 1:J_spec){
    y_spec[i] ~ dbinom(prob = ilogit(logit_spec[i]), size = n_spec[i]) 
  }
  # likelihood for 3 sens studies
  for(i in 1:J_sens){
    y_sens[i] ~ dbinom(prob = ilogit(logit_sens[i]), size = n_sens[i])
  }
  # likelihood for SC data
  p_sample <- p * ilogit(logit_sens[1]) + (1 - p) * (1 - ilogit(logit_spec[1])) # SC frequency of positive tests (prevalence corrected for imperfect tests)
  y_sample ~ dbinom(prob = p_sample, size = n_sample)
})

# inits <- list(
#   p = runif(n = 1),
#   mu_logit_spec = rnorm(n = 1,mean = 4,sd = 2),
#   mu_logit_sens = rnorm(n = 1,mean = 4,sd = 2),
#   sigma_logit_spec = rnorm(n = 1,mean = 0,sd = 1),
#   sigma_logit_sens = rnorm(n = 1,mean = 0,sd = 1),
#   logit_spec = rnorm(n = 14,mean = 4,sd = 2),
#   logit_sens = rnorm(n = 4,mean = 4,sd = 2)
# )

sc_hier_model_n <- nimble::nimbleModel(
  code = sc_hier_n,
  data = list(
    y_sample=50,
    y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
    y_sens=c(0, 78, 27, 25)
  ),
  constants = list(
    n_sample=3330,
    J_spec=14,
    n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
    J_sens=4,
    n_sens=c(0, 85, 37, 35)
  )
)

mon_vars <- c("p","mu_logit_spec","mu_logit_sens","sigma_logit_spec","sigma_logit_sens","logit_spec","logit_sens","p_sample")
sc_hier_model_mcmc_n <- nimble::buildMCMC(conf = sc_hier_model_n,monitors = mon_vars)

sc_hier_model_cpp_n <- nimble::compileNimble(sc_hier_model_n)

sc_hier_model_cpp_mcmc_n <- nimble::compileNimble(sc_hier_model_mcmc_n, project = sc_hier_model_cpp_n)

samples <- nimble::runMCMC(mcmc = sc_hier_model_cpp_mcmc_n,niter = 1e5,nburnin = 1e4,thin = 10,progressBar = TRUE)


hist(samples[,"p"],breaks = 50)
hist(extract(sc_hier_model,"p")[[1]],breaks = 50)


# sc_model_n <- nimble::nimbleModel(
#   code = santa_clara_n,
#   data = list(y_sample=50, y_spec=369+30, y_sens=25+78),
#   constants = list(n_sample = 3330, n_spec=401, n_sens=122),
#   inits = list(p = 0.1, spec = 0.9, sens= 0.9)
# )


# --------------------------------------------------------------------------------
#   even better model: hierarchical, MRP
# --------------------------------------------------------------------------------

sc_hier_mrp <- '
  data {
    int<lower = 0> N;  // number of tests in the sample (3330 for Santa Clara)
    int<lower = 0, upper = 1> y[N];  // 1 if positive, 0 if negative
    vector<lower = 0, upper = 1>[N] male;  // 0 if female, 1 if male
    int<lower = 1, upper = 4> eth[N];  // 1=white, 2=asian, 3=hispanic, 4=other
    int<lower = 1, upper = 4> age[N];  // 1=0-4, 2=5-18, 3=19-64, 4=65+
    int<lower = 0> N_zip;  // number of zip codes (58 in this case)
    int<lower = 1, upper = N_zip> zip[N];  // zip codes 1 through 58
    vector[N_zip] x_zip;  // predictors at the zip code level
    int<lower = 0> J_spec;
    int<lower = 0> y_spec [J_spec];
    int<lower = 0> n_spec [J_spec];
    int<lower = 0> J_sens;
    int<lower = 0> y_sens [J_sens];
    int<lower = 0> n_sens [J_sens];
    int<lower = 0> J;  // number of population cells, J = 2*4*4*58
    vector<lower = 0>[J] N_pop;  // population sizes for poststratification
    real<lower = 0> coef_prior_scale;
    real<lower = 0> logit_spec_prior_scale;
    real<lower = 0> logit_sens_prior_scale;
  }
  parameters {
    real mu_logit_spec;
    real mu_logit_sens;
    real<lower = 0> sigma_logit_spec;
    real<lower = 0> sigma_logit_sens;
    vector<offset = mu_logit_spec, multiplier = sigma_logit_spec>[J_spec] logit_spec;
    vector<offset = mu_logit_sens, multiplier = sigma_logit_sens>[J_sens] logit_sens;
    vector[3] b;  // intercept, coef for male, and coef for x_zip
    real<lower = 0> sigma_eth;
    real<lower = 0> sigma_age;
    real<lower = 0> sigma_zip;
    vector<multiplier = sigma_eth>[4] a_eth;  // varying intercepts for ethnicity
    vector<multiplier = sigma_age>[4] a_age;  // varying intercepts for age category
    vector<multiplier = sigma_zip>[N_zip] a_zip;  // varying intercepts for zip code
  }
  transformed parameters {
    vector[J_spec] spec = inv_logit(logit_spec);
    vector[J_sens] sens = inv_logit(logit_sens);
  }
  model {
    vector[N] p = inv_logit(b[1]
                            + b[2] * male
                            + b[3] * x_zip[zip]
                            + a_eth[eth]
                            + a_age[age]
                            + a_zip[zip]);
    vector[N] p_sample = p * sens[1] + (1 - p) * (1 - spec[1]);
    y ~ bernoulli(p_sample);
    y_spec ~ binomial(n_spec, spec);
    y_sens ~ binomial(n_sens, sens);
    logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
    logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
    sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
    sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
    mu_logit_spec ~ normal(4, 2);  // weak prior on mean of distribution of spec
    mu_logit_sens ~ normal(4, 2);  // weak prior on mean of distribution of sens
    a_eth ~ normal(0, sigma_eth);
    a_age ~ normal(0, sigma_age);
    a_zip ~ normal(0, sigma_zip);
    // prior on centered intercept
    b[1] + b[2] * mean(male) + b[3] * mean(x_zip[zip]) ~ logistic(0, 1);
    b[2] ~ normal(0, coef_prior_scale);
    b[3] ~ normal(0, coef_prior_scale / sd(x_zip[zip]));  // prior on scaled coef
    sigma_eth ~ normal(0, coef_prior_scale);
    sigma_age ~ normal(0, coef_prior_scale);
    sigma_zip ~ normal(0, coef_prior_scale);
  }
  generated quantities {
    real p_avg;
    vector[J] p_pop;  // population prevalence in the J poststratification cells
    int count;
    count = 1;
    for (i_zip in 1:N_zip) {
      for (i_age in 1:4) {
        for (i_eth in 1:4) {
          for (i_male in 0:1) {
            p_pop[count] = inv_logit(b[1]
                                     + b[2] * i_male
                                     + b[3] * x_zip[i_zip]
                                     + a_eth[i_eth]
                                     + a_age[i_age]
                                     + a_zip[i_zip]);
            count += 1;
          }
        }
      }
    }
    p_avg = sum(N_pop .* p_pop) / sum(N_pop);
  }
'

# To fit the model, we need individual-level data.  These data are not publcly available, so just to get the program running, we take the existing 50 positive tests and assign them at random to the 3330 people.
N <- 3330
y <- sample(rep(c(0, 1), c(3330 - 50, 50)))
n <- rep(1, 3330)

# Here are the counts of each sex, ethnicity, and age from Bendavid et al. (2020).  We don't have zip code distribution but we looked it up and there are 58 zip codes in Santa Clara County; for simplicity we asssume all zip codes are equally likely.  We then assign these traits to people at random.  This is wrong--actually, these variable are correlated in various ways--but, again,. now we have fake data we can use to fit the model.
male <- sample(rep(c(0,1), c(2101, 1229)))
eth <- sample(rep(1:4, c(2118, 623, 266, 306+17)))
age <- sample(rep(1:4, c(71, 550, 2542, 167)))
N_zip <- 58
zip <- sample(1:N_zip, 3330, replace=TRUE)

# Setting up the zip code level predictor.  In this case we will use a random number with mean 50 and standard deviation 20.  These are arbitrary numbers that we chose just to be able to test the centering and scaling in the model.   In real life we might use %Latino or average income in the zip code
x_zip <- rnorm(N_zip, 50, 20)

# Setting up the poststratification table.  For simplicity we assume there are 1000 people in each cell in the county.  Actually we'd want data from the Census.
J <- 2*4*4*N_zip
N_pop <- rep(NA, J)
count = 1
for (i_zip in 1:N_zip){
  for (i_age in 1:4){
    for (i_eth in 1:4){
      for (i_male in 0:1){
        N_pop[count] <- 1000
        count = count + 1
      }
    }
  }
}

# Put togther the data and fit the model
santaclara_mrp_data <- list(
  N=N,
  y=y,
  male=male,
  eth=eth,
  age=age,
  zip=zip,
  N_zip=N_zip,
  x_zip=x_zip,
  J_spec=14,
  y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
  n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
  J_sens=4,
  y_sens=c(0, 78, 27, 25),
  n_sens=c(0, 85, 37, 35),
  logit_spec_prior_scale=0.3,
  logit_sens_prior_scale=0.3,
  coef_prior_scale=0.5,
  J=J,
  N_pop=N_pop
)

sc_hier_mrp_model <- rstan::stan(
  model_code = sc_hier_mrp,
  data = santaclara_mrp_data,
  chains = 4,
  warmup = 1e4,
  iter = 2e4,
  cores = 4,
  control=list(adapt_delta=0.95),
  verbose = TRUE
)
