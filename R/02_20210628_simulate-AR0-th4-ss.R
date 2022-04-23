# load packages
library(tidyverse)
library(furrr)

# source functions
source(here::here("R/01a_20210507_simulateLongitudinalFunction.R"))
source(here::here("R/02a_20210511_twoStageMCEMFunction.R"))

# inverse logit function
logit = function(x) log(x/(1-x))

# set parameters
m=50; n=100; t=25; theta=4; sigma=0.5
formula.mu = list(x1 ~ 1, x2 ~ 1); formula.phi = list(~ 1, ~ 1); formula.p = list(~ 1, ~ 1)
covariates=NULL; AR = "AR0"; prec=200
rho = list(logit(0.4), logit(0.50)); delta = list(logit(2/7), logit(3/9)); kappa = list(log(7), log(9))
set.seed(758)

# simulate data
df.ar0 = lapply(1:m, function(i) {
  simulateLongitudinalZIBR(n, t, theta, sigma, formula.mu, formula.phi, formula.p, covariates, AR, rho, delta, kappa) 
})

# two-stage estimation
ts.mcem = vector("list", length=m)

for (i in seq_along(ts.mcem)) {
  ts.mcem[[i]] = twoStageMCEM(formula.mu=formula.mu, formula.phi=formula.phi, formula.p=formula.p, 
                              df=df.ar0[[i]][,-3], x=c("x1", "x2"), timepoints="tp", AR=AR, 
                              proposal.sd=rep(1,t), niter=4000, emruns=20, ncore=t)
  
  print(i)
}

# save data
saveRDS(df.ar0, file="20210916-simData-ar0-theta4-sigma-s.RDS")
saveRDS(ts.mcem, file="20210916-tsmcem-ar0-theta4-sigma-s.RDS")