
library(tmap)
library(tmaptools)
library(tigris)
library(tidyverse)
library(OpenStreetMap)
library(rJava)
library(osmdata)
library(cartography)
library(sp)
library(tidycensus)

library(nimble, warn.conflicts = FALSE)
library(CARBayesdata, quietly = TRUE)
library(sp, quietly = TRUE)
library(spdep, quietly = TRUE)

data(GGHB.IG)
data(respiratorydata)
#We handle the spatial analysis here with nb2WB from the package spdep.
respiratorydata_spatial <- merge(x = GGHB.IG, y = respiratorydata, 
                                 by = "IG", all.x = FALSE)
W.nb <- poly2nb(respiratorydata_spatial, row.names =  rownames(respiratorydata_spatial@data))
## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W.nb)

# A vector of indices indicating which regions are neighbors of which.
nbInfo$adj

stop()

# get zipcodes for covid study
setwd("~/Google Drive/02_research/covid/info_flyer1/")
db <- load("tmap_scrn_data.RData"); db
options(tigris_use_cache = TRUE)
zcta_bounds_region <- zctas(cb = T, starts_with = unique(map_scrn_zip$GEOID10))
zcta_bounds_region

W.nb <- poly2nb(zcta_bounds_region, 
                row.names = zcta_bounds_region$ZCTA5CE10)
nbInfo <- nb2WB(W.nb)

# A vector of indices indicating which regions are neighbors of which.
nbInfo$adj

# A vector of weights. In this case, all weights are 1.
head(nbInfo$weights)

# A vector of length N. num[n] indicates how many neighbors region n contains.
# This helps map the adj vector to the starting region.
length(nbInfo$num)
nbInfo$num


################################################################################  
# COVID Ex
################################################################################  

# set up data -------
N <- 3330
n <- rep(1, 3330)
#y <- sample(rep(c(0, 1), c(3330 - 50, 50)))

# Here are the counts of each sex, ethnicity, and age from Bendavid et al. (2020).  We don't have zip code distribution but we looked it up and there are 58 zip codes in Santa Clara County; for simplicity we asssume all zip codes are equally likely.  We then assign these traits to people at random.  This is wrong--actually, these variable are correlated in various ways--but, again,. now we have fake data we can use to fit the model.
male <- sample(rep(c(0,1), c(2101, 1229)))
eth <- sample(rep(1:4, c(2118, 623, 266, 306+17)))
age <- sample(rep(1:4, c(71, 550, 2542, 167)))
N_zip <- 31
zip <- sample(1:N_zip, 3330, replace=TRUE)

# sample y's dependent on zipco4
y <- rbinom(N, size = 1, prob = plogis((zip-mean(zip))/4)/24)
table(y, zip) %>% addmargins()
mean(y)

# Setting up the zip code level predictor.  In this case we will use a random number with mean 50 and standard deviation 20.  These are arbitrary numbers that we chose just to be able to test the centering and scaling in the model.   In real life we might use %Latino or average income in the zip code
x_zip <- rnorm(N_zip, 50, 20)

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


# NIMBLE version
sc_mrp_nimble <- nimble::nimbleCode({
  
  # prior models
  
  # hyperpriors for sd of random intercepts
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # ethnicity
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # age category
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # zip code
  
  # random intercepts
  for(i in 1:4){
    a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth)
    a_age[i] ~ dnorm(mean = 0, sd = sigma_age)
  }
  for(i in 1:N_zip){
    a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip)
  }
  
  # prior on centered intercept
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip_zip[]))
  # b[1] + b[2] * mean(male[]) + b[3] * mean(x_zip_zip[]) ~ dlogis(location = 0, scale = 1)
  b[1] ~ dlogis(location = 0 - b[2] * mean(male[]) - b[3] * mean(x_zip_zip[]), scale = 1)
  
  # likelihood model for SC data
  for(i in 1:N){
    p[i] <- ilogit(b[1] + b[2] * male[i] + b[3] * x_zip_zip[i] + 
                     a_eth[eth[i]] + a_age[age[i]] + a_zip[zip[i]]) 
    y[i] ~ dbern(prob = p[i])
  }
  
  # post-stratified prevalence (using https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays)

  for (i_zip in 1:N_zip) {
    for (i_age in 1:4) {
      for (i_eth in 1:4) {
        for (i_male in 0:1) {
          p_pop[(((i_zip-1) * 4*4*2) + ((i_age-1) * 4*2) + ((i_eth-1) * 2) + i_male) + 1] <- 
            ilogit(b[1] + b[2] * i_male + b[3] * x_zip[i_zip] +  a_eth[i_eth] + a_age[i_age] + a_zip[i_zip])
        }
      }
    }
  }
  p_avg <- inprod(N_pop[],p_pop[]) / sum(N_pop[])

})


sc_mrp_model_nimble <- nimble::nimbleModel(
  code = sc_mrp_nimble,
  dimensions = list(p_pop = J),
  data = list(y=y),
  constants = list(
    N=N,
    male=male,
    eth=eth,
    age=age,
    zip=zip,
    N_zip=N_zip,
    x_zip=x_zip,
    x_zip_zip=x_zip[zip],
    coef_prior_scale=0.5,
    J=J,
    N_pop=N_pop
  )
)
mon_vars <- c("p_avg", "p_pop")
sc_mrp_model_mcmc_n <- nimble::buildMCMC(conf = sc_mrp_model_nimble,monitors = mon_vars)
sc_mrp_model_nimble_cpp <- nimble::compileNimble(sc_mrp_model_nimble)
sc_mrp_model_mcmc_nimble_cpp <- nimble::compileNimble(sc_mrp_model_mcmc_n, 
                                                      project = sc_mrp_model_nimble_cpp)
sc_mrp_nsamp <- nimble::runMCMC(mcmc = sc_mrp_model_mcmc_nimble_cpp,niter = 1.1e5/4,nburnin = 1e4/4,thin = 10,progressBar = TRUE, summary = T)


################################
# dcar version
################################



# NIMBLE version
car_mrp_nimble <- nimble::nimbleCode({
  
  # prior models
  
  # car priors
  tau <- 1 / sigma^2
  
  # hyperpriors for sd of random intercepts
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # ethnicity
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # age category
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # zip code
  
  # random intercepts
  for(i in 1:4){
    a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth)
    a_age[i] ~ dnorm(mean = 0, sd = sigma_age)
  }
  for(i in 1:N_zip){
    a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip)
  }
  
  # prior on centered intercept
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip_zip[]))
  # b[1] + b[2] * mean(male[]) + b[3] * mean(x_zip_zip[]) ~ dlogis(location = 0, scale = 1)
  b[1] ~ dlogis(location = 0 - b[2] * mean(male[]) - b[3] * mean(x_zip_zip[]), scale = 1)
  
  # latent process
  s[1:N_zip] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N_zip], tau, zero_mean = 0)
  
  # likelihood model for SC data
  for(i in 1:N){
    p[i] <- ilogit(b[1] + b[2] * male[i] + b[3] * x_zip_zip[i] + 
                     a_eth[eth[i]] + a_age[age[i]] + a_zip[zip[i]] + s[zip[i]]) 
    y[i] ~ dbern(prob = p[i])
  }
  
  # post-stratified prevalence (using https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays)
  
  for (i_zip in 1:N_zip) {
    for (i_age in 1:4) {
      for (i_eth in 1:4) {
        for (i_male in 0:1) {
          p_pop[(((i_zip-1) * 4*4*2) + ((i_age-1) * 4*2) + ((i_eth-1) * 2) + i_male) + 1] <- 
            ilogit(b[1] + b[2] * i_male + b[3] * x_zip[i_zip] +  
                     a_eth[i_eth] + a_age[i_age] + a_zip[i_zip] + + s[i_zip])
        }
      }
    }
  }
  p_avg <- inprod(N_pop[],p_pop[]) / sum(N_pop[])
  
})

car_mrp_model_nimble <- nimble::nimbleModel(
  code = car_mrp_nimble,
  dimensions = list(p_pop = J),
  data = list(y=y),
  constants = list(
    N=N,
    male=male,
    eth=eth,
    age=age,
    zip=zip,
    N_zip=N_zip,
    x_zip=x_zip,
    x_zip_zip=x_zip[zip],
    coef_prior_scale=0.5,
    J=J,
    N_pop=N_pop, 
    # dcar data/constants
    L = length(nbInfo$adj), 
    adj = nbInfo$adj, 
    weights = nbInfo$weights, 
    num = nbInfo$num
  )
)


mon_vars <- c("p_avg", "p_pop", "tau")
car_mrp_model_mcmc_n <- nimble::buildMCMC(conf = car_mrp_model_nimble,monitors = mon_vars)
car_mrp_model_nimble_cpp <- nimble::compileNimble(car_mrp_model_nimble)
car_mrp_model_mcmc_nimble_cpp <- nimble::compileNimble(car_mrp_model_mcmc_n, 
                                                      project = car_mrp_model_nimble_cpp)
car_mrp_nsamp <- nimble::runMCMC(mcmc = sc_mrp_model_mcmc_nimble_cpp,niter = 1.1e5/4,nburnin = 1e4/4,thin = 10,progressBar = TRUE, summary = T)

sc_mrp_nsamp$summary %>% head
car_mrp_nsamp$summary %>% head

zip_index <- rep(1:31, J/31) 
zip_index %>% head
zip_index %>% tail
zip_index %>% length; J
zip_index

p_pop_car <- car_mrp_nsamp$samples[, -1]
p_pop_car[, zip_index==1] %>% rowMeans
p_pop_car_zip <- sapply(1:31, function(x) {
  p_pop_car[, zip_index==x] %>% rowMeans
})

p_pop_car_zip_summary <- apply(p_pop_car_zip, 2, function(x) {
  c(mean(x), quantile(x, probs = c(0.5, 0.025, 0.975)))
}) %>% data.frame() %>% t() %>% data.frame() %>% mutate(zip = 1:31) %>% select(zip, 1:4)

p_pop_mrp <- sc_mrp_nsamp$samples[, -1]
p_pop_mrp_zip <- sapply(1:31, function(x) {
  p_pop_mrp[, zip_index==x] %>% rowMeans
})
p_pop_mrp_zip_summary <- apply(p_pop_mrp_zip, 2, function(x) {
  c(mean(x), quantile(x, probs = c(0.5, 0.025, 0.975)))
}) %>% data.frame() %>% t() %>% data.frame() %>% mutate(zip = 1:31) %>% select(zip, 1:4)


summary <- full_join(p_pop_car_zip_summary, p_pop_mrp_zip_summary, by = "zip", suffix = c("mrp", "car"))

summary %>% head

# truth
table(zip, y) %>% prop.table(1) %>% data.frame() %>% filter(y==1) %>% 
  mutate(zip = as.integer(as.character(zip))) %>% 
  full_join(summary)


summary %>% ggplot() + 
  geom_point(aes(x=V1mrp, y = V1car)) + 
  geom_abline(slope=1)
