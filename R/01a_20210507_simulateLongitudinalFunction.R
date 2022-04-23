simulateLongitudinalZIBR = 
  function(n, t, theta, sigma, formula.mu, formula.phi, formula.p, covariates=NULL, AR = "AR0", rho, delta, kappa, prec=200){
    theta.t = rnorm(t, theta, sigma)
    
    source(here::here("R/01a_20210304_simulateDataFunction.R"))
    
    if(AR == "AR0"){
      simX.ls = lapply(theta.t, simulateBivariateZIBR, n=n, formula.mu=formula.mu, 
                       formula.phi=formula.phi, formula.p=formula.p, covariates=covariates, 
                       rho=rho, delta=delta, kappa=kappa, prec=prec)
    } else if(AR == "AR1"){
      simX.ls = vector(mode = "list", length = t)
      
      simX.ls[[1]] = simulateBivariateZIBR(n=n, theta=theta.t[[1]], formula.mu=formula.mu, 
                                           formula.phi=formula.phi, formula.p=formula.p, covariates=covariates, 
                                           rho=rho,delta=delta, kappa=kappa, prec=prec)
      
      for (s in 2:t) {
        covariates = simX.ls[[s-1]]
        colnames(covariates) = c("xtm1.1", "xtm1.2")
        
        formula.p = list(~ xtm1.1, ~ xtm1.2)
        formula.mu = list(x1 ~ xtm1.1, x2 ~ xtm1.2)
        formula.phi = list(~ 1, ~ 1)
        
        rho = list(c(1, -4), c(1, -4))
        delta = list(c(-1, 1.5), c(-1, 1.5))
        kappa = list(c(1), c(1))
        
        simX.ls[[s]] = simulateBivariateZIBR(n=n, theta=theta.t[[s]], formula.mu=formula.mu, 
                                             formula.phi=formula.phi, formula.p=formula.p, covariates=covariates, 
                                             rho = rho,delta =delta, kappa=kappa, prec=prec)
      }
      
    } else{
      stop("ERROR: AR must be equal to 'AR0' or 'AR1'.")
    }
    
    out.ls = lapply(1:t, function(tt) {data.frame(id = rownames(simX.ls[[tt]]), 
                                                  tp = rep(tt, n), 
                                                  theta.t = rep(theta.t[[tt]], n), 
                                                  simX.ls[[tt]])} )
    
    out = do.call(rbind, out.ls)
    
    return(out)
  }