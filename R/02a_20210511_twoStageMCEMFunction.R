# Frank Copula Theta to Kendall's Tau
frankKendallTau <- function(theta){
  1 - (4/theta) + (4/(theta^2))*(integrate(f = function(t){t/(exp(t) - 1)}, lower = 0, upper = theta)[[1]])
}

# Kendall's Tau to Frank Copula Theta
# Adapted from the VineCopula package (https://cran.r-project.org/web/packages/VineCopula/index.html)
# https://github.com/tnagler/VineCopula/blob/master/R/BiCopTau2Par.R
frankTau2Theta <- function(tau) {
  if (abs(tau) > 0.99999) {
    return(Inf)
  }
  
  if(tau > 0){
    a = 1
  } else if(tau < 0) {
    a = -1
    tau = -tau
  } else if(tau == 0){
    return(0)
  }
  
  v = uniroot(function(x) tau - frankKendallTau(x), 
              lower = 0 + .Machine$double.eps^0.5, upper = 100, 
              tol = .Machine$double.eps^0.5)$root
  
  return(a*v)
}

# Marginal MLE function with AR model
marginalMLE = function(formula.mu, formula.phi = ~ 1, formula.p = ~ 1, df, x, timepoints = NULL, AR){
  
  marginalMLE.t = function(formula.mu, formula.phi = ~ 1, formula.p = ~ 1, df){
    # make inputted formula of class "formula"
    formula.mu = as.formula(formula.mu); formula.phi = as.formula(formula.phi); formula.p = as.formula(formula.p)
    
    # fit ZI-Beta regression model
    model = suppressWarnings( # suppress function warnings
      gamlss::gamlss(formula = formula.mu, sigma.formula = formula.phi, nu.formula = formula.p, 
                     family = gamlss.dist::BEZI, data = df, 
                     control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
    )
    
    # check convergence
    convergenceFail = model$converged==FALSE 
    
    # if model failed to converge refit
    if(convergenceFail){
      
      message("Model failed to converge. Refitting.")
      
      model = gamlss::refit(model)
      
      if(model$converged==FALSE){
        message("Model refit and failed to converged.")
      }
    }
    
    # combine estimated regression coefficients 
    zibrCoef = unlist(gamlss::coefAll(model))
    names(zibrCoef) = janitor::make_clean_names(names(zibrCoef))
    
    # estimated P(x=1), mean, and dispersion
    p = fitted(model,"nu"); mu = fitted(model,"mu"); phi = fitted(model,"sigma")
    
    # estimated shape parameters
    alpha = mu*phi; beta = phi - mu*phi
    
    # output
    out = list("coefficients" = zibrCoef, 
               "fitted" = data.frame(p=p, mu=mu, phi=phi), 
               "shape" = data.frame(p=p, alpha=alpha, beta=beta))
    
    return(out)
  }
  
  if(is.null(timepoints)){ # if cross-sectional data run single marginalMLE
    out = marginalMLE.t(formula.mu, formula.phi, formula.p, df)
  } else{ # if longitudinal data repeat for each timepoint
    
    if(AR == "AR0"){ # if ar0 model estimate margins directly
      # split data by timepoint
      dfT.ls = split(df, df[,timepoints])
      
      # MLE/timepoint
      MLE.ls = lapply(dfT.ls, marginalMLE.t, formula.mu=formula.mu, formula.phi=formula.phi, formula.p=formula.p)
    } else if(AR == "AR1"){ # if ar1, format data to estimate margins
      # non time-varying variables
      idvars = colnames(df)[-which(colnames(df) %in% c(timepoints,x))]
      
      # long to wide data set
      df.wide = reshape(df, idvar = idvars, timevar = timepoints, direction = "wide")
      
      # unique values of the time indicator variable
      tp.id = unique(df[,timepoints])
      
      # split data by timepoint
      dfT.ls = vector("list", length(tp.id))                        # empty list for dfs
      names(dfT.ls) = unique(tp.id)
      dfT.ls[[1]] = df.wide[,c(idvars, paste0(x, ".", tp.id[[1]]))] # first time point no AR
      colnames(dfT.ls[[1]]) = c(idvars, x)                          # rename variables
      
      dfT.ls[names(dfT.ls)[-1]] = lapply(names(dfT.ls)[-1], function(i) {
        # covars = paste0(rep(x,2), ".", c(i, i-1))    # AR variables
        current.t = names(dfT.ls)[which(names(dfT.ls) == i)]
        previous.t = names(dfT.ls)[which(names(dfT.ls) == i)-1]
        covars = paste0(rep(x,2), ".", c(current.t, previous.t)) # AR variables
        dfT.ls[[i]] = df.wide[,c(idvars,covars)]                 # select all variables
        colnames(dfT.ls[[i]]) = c(idvars, x, "xtm1")             # rename
        return(dfT.ls[[i]])
      })
      
      # MLE/timepoint
      MLE.ls = vector("list", length(tp.id)) # empty list for dfs
      names(MLE.ls) = 1:length(MLE.ls)
      MLE.ls[[1]] = marginalMLE.t(formula.mu=formula.mu, formula.phi=formula.phi, formula.p=formula.p, df=dfT.ls[[1]])
      
      # update formulas with AR
      formula.mu  = update.formula(formula.mu, ~ . + xtm1)
      formula.phi = update.formula(formula.phi, ~ . + xtm1)
      formula.p   = update.formula(formula.p, ~ . + xtm1)
      
      MLE.ls[2:length(tp.id)] = lapply(dfT.ls[2:length(tp.id)], marginalMLE.t, 
                                       formula.mu=formula.mu, formula.phi=formula.phi, formula.p=formula.p)
      
      # add NA for x1
      mu.id = grep("mu", names(MLE.ls[[1]]$coefficients))
      phi.id = grep("phi", names(MLE.ls[[1]]$coefficients))
      nu.id = grep("nu", names(MLE.ls[[1]]$coefficients))
      
      MLE.ls[[1]]$coefficients = c(MLE.ls[[1]]$coefficients[mu.id], "mu_xtm1" = NA, 
                                   MLE.ls[[1]]$coefficients[phi.id], "phi_xtm1" = NA, 
                                   MLE.ls[[1]]$coefficients[nu.id], "nu_xtm1" = NA )
    } else{
      stop("ERROR: AR must be equal to 'AR0' or 'AR1'.")
    }
    
    # bind each timepoint MLE list into one list. Code from:
    # https://stackoverflow.com/questions/41596326/rbind-dataframes-across-nested-lists
    out = do.call(Map, c(f = rbind, MLE.ls))
    
    # add time variable
    out[2:3] = lapply(out[2:3], function(x){
      df = cbind(gsub("\\..*","", rownames(x)), x) # add timepoint var from rownames (before ".")
      rownames(df) = NULL                          # remove rownames
      colnames(df)[1] = timepoints                 # new var = timepoint
      return(df)
    })
  }
  
  return(out)
}

# Metropolis-Hastings Algorithm to sample \theta^{(t)}
thetaMH = function(theta, sigma, shape, x.df, timepoints, init=NULL, proposal.sd, niter=5000, ncore=1){
  
  tp = factor(shape[,timepoints], levels = unique(shape[,timepoints]))
  xT.ls = split(x.df, tp)
  shapeT.ls = split(shape, tp)
  
  t = length(shapeT.ls)
  
  theta.t.MH = function(theta, sigma, shape, x.df, timepoints, init, proposal.sd, niter){
    
    theta.t = vector(mode = "numeric", length = niter)
    accept = vector(mode = "numeric", length = niter-1)
    
    if(is.null(init)){
      theta.t[[1]] = theta
    } else{
      theta.t[[1]] = init
    }
    
    # observed data
    x1 = x.df[,1]; x2 = x.df[,2]
    # marginal shape parameters
    p1 = shape[["p1"]]; alpha1 = shape[["alpha1"]]; beta1 = shape[["beta1"]]
    p2 = shape[["p2"]]; alpha2 = shape[["alpha2"]]; beta2 = shape[["beta2"]]
    
    # source the distribution functions
    source(here::here("R/00b_20200304_distributionFunctions.R"), local = TRUE)
    
    # density under initial (current) value
    if(theta.t[[1]]==0){
      f.current = prod(fx(x1, p1, alpha1, beta1)*fx(x2, p2, alpha2, beta2))*dnorm(theta.t[[1]], mean = theta, sd = sigma)
    } else{
      f.current = prod(fxy(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, theta.t[[1]]))*dnorm(theta.t[[1]], mean = theta, sd = sigma)
    }
    
    for (i in 2:niter) {
      
      # draw candidate theta.t
      candidate = rnorm(1, mean = theta.t[[i-1]], sd = proposal.sd)
      
      # density under candidate value
      if(candidate==0){
        f.candidate = prod(fx(x1, p1, alpha1, beta1)*fx(x2, p2, alpha2, beta2))*dnorm(candidate, mean = theta, sd = sigma)
      } else {
        f.candidate = prod(fxy(x1, p1, alpha1, beta1, x2, p2, alpha2, beta2, candidate))*dnorm(candidate, mean = theta, sd = sigma)
      }
      
      # ratio
      r = min(f.candidate/f.current, 1)
      
      u = runif(1)
      
      if(u < r){ # accept with probability r
        theta.t[i] = candidate    # new value = candidate
        f.current = f.candidate   # f.candidate become f.current
        accept[i]=1
      } else{ # reject with probability 1-r
        theta.t[i] = theta.t[i-1] # otherwise no change: new value = current value
      }
    }
    
    print(paste0("Acceptance Rate: ", sum(accept)/length(accept)*100, "%"))
    
    return(theta.t)
  }
  
  theta.t = vector("list", length = t)
  
  if(ncore==1){
    for (i in 1:t) {
      theta.t[[i]] = theta.t.MH(theta=theta, sigma=sigma, shape=shapeT.ls[[i]], x.df=xT.ls[[i]], 
                                init=init, proposal.sd=proposal.sd[[t]], niter=niter)
    }
  } else{
    ncore = min(ncore, future::availableCores()-1)
    future::plan(future::multisession, workers = ncore)
    
    theta.t = furrr::future_pmap(list(shapeT.ls,xT.ls,proposal.sd), function(a,b,c) {
      theta.t.MH(theta=theta, sigma=sigma, shape=a, x.df=b, init=init, proposal.sd=c, niter=niter)
      }, .options = furrr::future_options(seed = T) )
  }
  
  return(theta.t)
}

# Monte-Carlo EM Algorithm
MCEM = function(shape, x.df, timepoints, init=NULL, proposal.sd, niter=5000, emruns=25, ncore=1){
  theta.t = vector("list", length=emruns)
  theta = vector("list", length=emruns)
  sigma.sq = vector("list", length=emruns)
  sigma = vector("list", length=emruns)
  
  t = length(unique(shape[,timepoints]))
  
  for(i in 1:emruns) {
    
    # MCMC
    if(i==1){
      # theta0 = rnorm(1)
      tau0 = cor(x.df[,1], x.df[,2], method="kendall")
      theta0 = frankTau2Theta(tau0)
      sigma0 = abs(rnorm(1, 2))
      
      theta.t[[i]] = thetaMH(theta=theta0, sigma=sigma0, shape=shape, x.df=x.df, timepoints = timepoints, 
                             init=init, proposal.sd=proposal.sd, niter=niter, ncore)
    } else{
      theta.t[[i]] = thetaMH(theta=theta[[i-1]], sigma=sigma[[i-1]], shape=shape, x.df=x.df,
                             timepoints=timepoints, init=init, proposal.sd=proposal.sd, niter=niter, ncore)
    }
    
    # EM -- theta
    theta.t.bar = unlist(lapply(theta.t[[i]], mean))
    theta[[i]] = mean(theta.t.bar)
    
    # EM -- sigma
    theta.t.center.sq.bar = unlist(lapply(theta.t[[i]], function(k) mean((k - theta[[i]])^2)))
    sigma.sq[[i]] = (1/t)*sum(theta.t.center.sq.bar)
    sigma[[i]] = sqrt(sigma.sq[[i]])
    
    proposal.sd = 2.38*unlist(lapply(theta.t[[i]], sd))
    # proposal.sd = sd(theta.t.bar)
  }
  
  theta = unlist(theta)
  sigma.sq = unlist(sigma.sq)
  sigma = unlist(sigma)
  
  return(list(theta.t = theta.t, 
              theta = theta, 
              sigma = sigma #, 
              # prop.sd = proposal.sd
              ))
  
}

twoStageMCEM = function(formula.mu, formula.phi, formula.p, df, x, timepoints=NULL, AR, init=NULL,proposal.sd, niter, emruns, ncore){
  # marginal parameters
  marginalMLEs1 = marginalMLE(formula.mu[[1]], formula.phi[[1]], formula.p[[1]], df[,-which(names(df) == x[2])], x[1], timepoints, AR)
  marginalMLEs2 = marginalMLE(formula.mu[[2]], formula.phi[[2]], formula.p[[2]], df[,-which(names(df) == x[1])], x[2], timepoints, AR)
  
  # fix names
  # names(marginalMLEs1$coefficients) = paste0(names(marginalMLEs1$coefficients), 1)
  names(marginalMLEs1$fitted) = paste0(names(marginalMLEs1$fitted), 1)
  names(marginalMLEs1$shape) = paste0(names(marginalMLEs1$shape), 1)
  
  # names(marginalMLEs2$coefficients) = paste0(names(marginalMLEs2$coefficients), 2)
  names(marginalMLEs2$fitted) = paste0(names(marginalMLEs2$fitted), 2)
  names(marginalMLEs2$shape) = paste0(names(marginalMLEs2$shape), 2)
  
  names(marginalMLEs1) = paste0(names(marginalMLEs1), 1)
  names(marginalMLEs2) = paste0(names(marginalMLEs2), 2)
  
  shape = cbind(marginalMLEs1$shape1, marginalMLEs2$shape2)
  shape = shape[,-5]
  colnames(shape)[1] = paste0(timepoints)
  
  x.df = df[,x]
  
  # MCEM
  theta.mcem = MCEM(shape = shape, x.df = x.df, timepoints = timepoints, init, proposal.sd, niter, emruns, ncore)
  
  return(list(margins1=marginalMLEs1, 
              margins2=marginalMLEs2, 
              mcem=theta.mcem))
  }

