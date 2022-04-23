# load packages
library(tidyverse)
library(furrr)

# set up parallelization
nCores = 50
future::plan(future::multisession, workers = nCores)

# source functions
source(here::here("R/01a_20210507_simulateLongitudinalFunction.R"))
source(here::here("R/02a_20210511_twoStageMCEMFunction.R"))
source(here::here("R/03a_20210726_MCLRTFunction.R"))

# inverse logit function
logit = function(x) log(x/(1-x))

# set parameters
m=500; n=100; t=25; theta=0; sigma=0.5
formula.mu = list(x1 ~ 1, x2 ~ 1); formula.phi = list(~ 1, ~ 1); formula.p = list(~ 1, ~ 1)
covariates=NULL; AR = "AR0"; prec=200
rho = list(logit(0.4), logit(0.50)); delta = list(logit(2/7), logit(3/9)); kappa = list(log(7), log(9))
set.seed(946)

# simulate data
df.ar0 = lapply(1:m, function(i) {
  simulateLongitudinalZIBR(n, t, theta, sigma, formula.mu, formula.phi, formula.p, covariates, AR, rho, delta, kappa) 
  })

# two-stage estimation
df = lapply(df.ar0, function(x) x[,-3])

ts.mcem = furrr::future_map(df, twoStageMCEM, 
                            formula.mu=formula.mu, formula.phi=formula.phi, formula.p=formula.p, 
                            x=c("x1", "x2"), timepoints="tp", AR=AR, proposal.sd=rep(1,t), 
                            niter=2500, emruns=12, ncore=1, 
                            .progress = T, .options = furrr::future_options(seed = T))

out = tibble::tibble(df = df.ar0, 
                     tsMCEM = ts.mcem) %>% 
  dplyr::mutate(., x.df = furrr::future_map(df, function(df, x){df[,x]}, x=c("x1", "x2")) ) %>% 
  tidyr::unnest_wider(tsMCEM) %>% tidyr::unnest_wider(margins1) %>% tidyr::unnest_wider(margins2) %>% # unnest
  dplyr::mutate(., shape = purrr::map2(shape1, shape2, cbind)) %>% # single shape variable
  dplyr::mutate(., shape = purrr::map(shape, function(df) {        # edit shape var colnames
    df = df[,-which(colnames(df)=="tp2")]                          # drop 2nd timepoint (month) var
    colnames(df)[which(colnames(df)=="tp1")] = "tp"                # rename timepoint (month) var
    return(df) }) ) %>% 
  dplyr::mutate(., thetaNull = 0) %>% # null value
  dplyr::mutate(., thetaMLE = purrr::map_dbl(mcem, function(x) mean(tail(x$theta,5)) ), # theta MLE
                sigmaMLE = purrr::map_dbl(mcem, function(x) mean(tail(x$sigma,5)) ) )   # sigma MLE

# mcLRT
out = out %>% 
  dplyr::mutate(., mc.LRT = furrr::future_pmap(dplyr::select(., thetaNull,thetaMLE,sigmaMLE,shape,df), 
                                               mcLRT, 
                                               x=c("x1", "x2"), timepoints="tp", init=0, 
                                               proposal.sd=rep(1,t), niter=1500, emruns=12, ncore=1, 
                                               .progress = T, .options = furrr::future_options(seed = T)) )

# save data
saveRDS(out, file="20211217_simData-ar0-theta0-sigma-s.RDS")
