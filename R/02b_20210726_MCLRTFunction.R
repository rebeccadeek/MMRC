# Monte Carlo Likelihood Ratio Test
mcLRT = function(thetaNull, thetaMLE, sigmaMLE, shape, df, x, timepoints, init=NULL, proposal.sd, niter=5000, emruns=25, ncore=1){
  # Calculate sigmaNull
  theta.t  = vector("list", length=emruns)
  sigma.sq = vector("list", length=emruns)
  sigma    = vector("list", length=emruns)
  
  t = length(unique(shape[,timepoints]))
  
  x.df = df[,x]
  
  # source the distribution functions
  source(here::here("code/02a_20210511_twoStageMCEMFunction.R"), local = TRUE)
  
  for(i in 1:emruns) {
    
    # MCMC
    if(i==1){
      sigma.init = abs(rnorm(1, 2))
      
      theta.t[[i]] = thetaMH(theta=thetaNull, sigma=sigma.init, shape=shape, x.df=x.df, timepoints = timepoints, 
                             init=init, proposal.sd=proposal.sd, niter=niter, ncore)
    } else{
      theta.t[[i]] = thetaMH(theta=thetaNull, sigma=sigma[[i-1]], shape=shape, x.df=x.df,
                             timepoints=timepoints, init=init, proposal.sd=proposal.sd, niter=niter, ncore)
    }
    
    # EM -- sigma
    theta.t.center.sq.bar = unlist(lapply(theta.t[[i]], function(k) mean((k - thetaNull)^2)))
    sigma.sq[[i]] = (1/t)*sum(theta.t.center.sq.bar)
    sigma[[i]] = sqrt(sigma.sq[[i]])
    
    proposal.sd = 2.38*unlist(lapply(theta.t[[i]], sd))
  }
  
  sigma = unlist(sigma)
  sigmaNull = mean(tail(sigma, 5))
  
  # MH sampling of \theta^{(t)} under the null
  theta.t.null = thetaMH(theta=thetaNull, sigma=sigmaNull, shape=shape, x.df=x.df, timepoints=timepoints, 
                         init=init, proposal.sd=proposal.sd, niter=niter, ncore=ncore)
  
  theta.t.null.unls = as.numeric(unlist(theta.t.null))
  
  # Monte-Carlo Likelihood Ratio
  # source the distribution functions
  source(here::here("code/00b_20200304_distributionFunctions.R"), local = TRUE)
  
  # split data
  tp = factor(shape[,timepoints], levels = unique(shape[,timepoints]))
  xT.ls = split(x.df, tp)
  shapeT.ls = split(shape, tp)
  
  fxy_ifelse = function(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta){
    if(theta[[1]]==0){
      fx(x1, p1, alpha1, beta1)*fx(x2, p2, alpha2, beta2)
    } else {
      fxy(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta)
    }
  }
  ###################################################################################################################
  f.theta.mle  = lapply(theta.t.null, function(x) sum(dnorm(x, mean = thetaMLE, sd = sigmaMLE)))
  log.f.theta.mle = sum(log(unlist(f.theta.mle) ) )
  
  f.theta.null = lapply(theta.t.null, function(x) sum(dnorm(x, mean = thetaNull, sd = sigmaNull)))
  log.f.theta.null = sum(log(unlist(f.theta.null) ) )
  
  LR.ts = 2*(log.f.theta.mle - log.f.theta.null)
  
  return(list(sigma = sigma, 
              sigmaNull = sigmaNull, 
              theta.t.null = theta.t.null, 
              LR.ts = LR.ts))
}