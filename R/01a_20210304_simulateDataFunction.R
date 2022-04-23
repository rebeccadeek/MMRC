# Function to simulate correlated bivariate data using the Frank copula
# When theta != 0: simulate u and w from independent U(0,1) distributions where 
# u is the CDF for x1 and w is the conditional Frank copula for u and v (the CDF of x2)
# use the Rosenblatt transformation to get find v. Then use the probability integral
# transformation (PIT) to find x1 and x2
# When theta = 0 simulate u and v from independent U(0,1) and use the PIT to find x1 and x2
simulateBivariateZIBR = function(n, theta, formula.mu, formula.phi, formula.p, covariates=NULL, rho, delta, kappa, prec = 200) {
  # 0. convert covariates to data.matrix so can do %*%
  if(!is.null(covariates)) covariates = data.matrix(covariates)
  
  # 1. Transformation from rho/delta/kappa to p_i/alpha_i/beta_i
  
  # assign regression parameters their own vectors
  rho1 = rho[[1]]; rho2 = rho[[2]]
  delta1 = delta[[1]]; delta2 = delta[[2]]
  kappa1 = kappa[[1]]; kappa2 = kappa[[2]]
  
  # assign separate design matrices from covariate ls() & add column of 1s for intercept
  # Q (p)
  if(formula.p[[1]][[2]]==1){ # if RHS of y ~ x is == ~ 1 then intercept only model
    Q1 = cbind(rep(1, n))
  } else { # else extract covariates in the model
    Q1 = cbind("int"=rep(1, n), covariates[,all.vars(formula.p[[1]][[2]])])
  }
  
  if(formula.p[[2]][[2]]==1){
    Q2 = cbind(rep(1, n))
  } else {
    Q2 = cbind("int"=rep(1, n), covariates[,all.vars(formula.p[[2]][[2]])])
  }
  
  # W (mu)
  if(formula.mu[[1]][[3]]==1){
    W1 = cbind(rep(1, n))
  } else {
    W1 = cbind("int"=rep(1, n), covariates[,all.vars(formula.mu[[1]][[3]])])
  }
  
  if(formula.mu[[2]][[3]]==1){
    W2 = cbind(rep(1, n))
  } else {
    W2 = cbind("int"=rep(1, n), covariates[,all.vars(formula.mu[[2]][[3]])])
  }
  
  # Z (phi)
  if(formula.phi[[1]][[2]]==1){
    Z1 = cbind(rep(1, n))
  } else {
    Z1 = cbind("int"=rep(1, n), covariates[,all.vars(formula.phi[[1]][[2]])])
  }
  
  if(formula.phi[[2]][[2]]==1){
    Z2 = cbind(rep(1, n))
  } else {
    Z2 = cbind("int"=rep(1, n), covariates[,all.vars(formula.phi[[2]][[2]])])
  }
  
  # ZI-Beta Regression
  # GLM for p
  logit.p1 = Q1 %*% rho1; logit.p2 = Q2 %*% rho2
  # GLM for mu
  logit.mu1 = W1 %*% delta1; logit.mu2 = W2 %*% delta2
  # GLM for phi
  log.phi1 = Z1 %*% kappa1; log.phi2 = Z2 %*% kappa2
  
  # inverse logit function
  inv.logit = function(x) exp(x) / (1 + exp(x))
  
  # p
  p1 = inv.logit(logit.p1); p2 = inv.logit(logit.p2)
  # mu
  mu1 = inv.logit(logit.mu1); mu2 = inv.logit(logit.mu2)
  # phi
  phi1 = exp(log.phi1); phi2 = exp(log.phi2)
  
  # shape parameters
  alpha1 <- mu1*phi1; beta1 <- phi1 - mu1*phi1
  alpha2 <- mu2*phi2; beta2 <- phi2 - mu2*phi2
  
  # 2. Simulation bivariate data ZIB data
  
  # Checks
  ## initialize check for x1 < 3 & x2 < 3 as TRUE
  # check0 = TRUE
  checkX0 = FALSE; checkNZPairs = FALSE
  ## check if theta == 0 (independence)
  thetaCheck = theta != 0
  
  # simulate x1 and x2 using ZIB regression margins
  while(!(checkX0 & checkNZPairs)) {
    # simulate u and w, independent uniform random variables
    u = runif(n = n, min = 0, max = 1)                                      # u = F1(x1)
    
    # if theta != 0 simulate v by:
    if(thetaCheck){
      w = runif(n = n, min = 0, max = 1)                                    # w = c1(u,v)
      
      a = Rmpfr::mpfr(exp(-theta), prec)                                    # parameterize a = exp(-theta)
      v = as.numeric((-1/theta) * log(1 + (w*(a - 1))/(w + (a^u)*(1 - w)))) # v = F2(x2)
    }
    
    # if theta == 0 simulate v by:
    else {
      v = runif(n = n, min = 0, max = 1)                                    # v = F2(x2)
    }
    
    # **the inverse CDF transformation is done in two stages because doing it in one qbeta gives a warning**
    # if u <= p1, x1u = 0, else transform u into x1u using (u - p1)/(1 - p1)
    # likewise for v and x2v
    x1u = dplyr::if_else(u <= p1, 0, (u-p1)/(1-p1))
    x2v = dplyr::if_else(v <= p2, 0, (v-p2)/(1-p2))
    
    # if x1u == 0, x1 = 0, else x1 is the inverse CDF ox x1u
    # likewise for x2v and x2
    x1 = dplyr::if_else(x1u == 0, 0, qbeta(x1u, alpha1, beta1))      # use the inverse CDF to solve for x1
    x2 = dplyr::if_else(x2v == 0, 0, qbeta(x2v, alpha2, beta2))      # use the inverse CDF to solve for x2
    
    # check0 = length(which(x1 != 0)) < 3 | length(which(x2 != 0)) < 3 # each x needs >= 3 non-zero obs to estimate parameters
    checkX0 = length(which(x1 != 0)) > 3 | length(which(x2 != 0)) > 3 # each x needs >= 3 non-zero obs to estimate parameters
    checkNZPairs = length(which(x1!=0 & x2!=0)) > 1                   # at least 2 non-zero pairs are needed to stable est
  }
  
  return(x = data.frame(x1, x2))
}

# Function for df of different parameter settings for simulation
paramsFunction = function(n, theta) {
  # inverse logit function
  logit = function(x) log(x/(1-x))
  
  rho = list(
    list(logit(0.1), logit(0.25) ), 
    list(logit(0.4), logit(0.50) ), 
    list(logit(0.6), logit(0.75) ), 
    list(logit(0.2), logit(0.70) )
  )
  delta = list(
    list(logit(2/7), logit(3/9) ), 
    list(logit(1/2), logit(1/3) ), 
    list(logit(5/7), logit(2/3) ), 
    list(logit(2/7), logit(2/3) ) 
  )
  kappa = list(
    list(log(7), log(9) ), 
    list(log(4), log(6) ), 
    list(log(7), log(9) ), 
    list(log(7), log(9) )
  )
  
  param = tibble::tibble(n = rep(n, times = 4*4), 
                         theta = rep(theta, times = 4*4),
                         rho = rep(rho, times = 4), 
                         delta = rep(delta, each = 4), 
                         kappa = rep(kappa, each = 4)
  )
  
  return(param)
}

paramsFunctionQ1 = function(n, theta) {
  param = tibble::tibble(n = n, theta = theta,
                         rho = list(
                           list(c(-0.5, 0.7), c(-0.3, 0.4) ), 
                           list(c(-0.1, 0.7), c(0.1, 0.4) ),
                           list(c(0.5, 0.7), c(0.8, 0.4) ) 
                           ), 
                         delta = list(-0.7, -1) %>% list(),
                         kappa = list(1.5, 1.5) %>% list()
                         )
  
  return(param)
}
