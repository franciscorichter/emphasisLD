speciation_rate <- function(tm,
                            tree,
                            diversification_model,
                            sum_rates=FALSE){
  speciation_r = get(paste0("lambda.", diversification_model$model))
  lambda = speciation_r(tm,tree,diversification_model$pars,sum_rates=sum_rates)
  return(lambda)
}

extinction_rate <- function(tm,
                            tree,
                            diversification_model=NULL,
                            sum_rates=FALSE){
  mu = mu.constant(tm = tm,
                   tree = tree,
                   pars = diversification_model$pars,
                   sum_rates = sum_rates)
  return(mu)
}

sum_of_rates <- function(tm,
                         tree,
                         diversification_model){
  val = speciation_rate(tm,
                  tree,
                  diversification_model,
                  TRUE) + 
        extinction_rate(tm,
                   tree,
                   diversification_model,
                   TRUE)
  return(val)
}

##########################################################
##### lineages dependent speciation rate #################
##########################################################

lambda.ldpd <- function(tm,
                        tree,
                        pars,
                        sum_rates=FALSE){
  tree_extant = get_extant(tm,tree)
  gpd = GPD(tm,tree_extant)
  N = nrow(gpd)
  lambdas = rep(max(0,pars[2] + pars[3]*N),times=N) + pars[4]*colSums(gpd)/(N-1)
  if(sum_rates) lambdas = sum(lambdas)
  return(lambdas)
}

mu.constant <- function(tm,
                        tree,
                        pars,
                        sum_rates=FALSE){
  mu = max(0,pars[1])
  N = sapply(tm, n_from_time,tree=tree)
  if(sum_rates){
    mu = N*mu
  }else{
    mu = rep(mu,N)
  }
  return(mu)
}


nh_speciation_rate <- function(tm,tree,diversification_model){
  
  speciation_rate(tm,
                  tree,
                  diversification_model,
                  sum_rates=TRUE)*(1-exp(-diversification_model$pars[1]*(tree$ct-tm)))

}

# Speciations rates 

lambda.rpd1 <- function(tm,tree,pars,sum_rates=FALSE){
  N = sapply(tm, n_from_time,tree=tree)
  lambda = rep(max(0, pars[2] + pars[3]*N ),N)
  if(sum_rates) lambda <- sum(lambda)
  return(lambda)
}

lambda.rpd5c <- function(tm,tree,pars,sum_lambda=FALSE){
  
  pd = sapply(tm,phylodiversity,tree=tree)-tm
  N = sapply(tm, n_from_time,tree=tree)
  lambda = rep(max(0, pars[2] + pars[3]*N + pars[4] * pd/N ) ,N)
  if(sum_lambda) lambda <- sum(lambda)
  return(lambda)
  
}

nh_rate <- function(x,model,pars,tree){
  speciation_rate(tm=x,
                  tree = tree,
                  pars = pars,
                  model = model,
                  sum_lambda = TRUE)*
    (1-exp( -(tree$ct-x) * extinction_rate(pars = pars) ))
}
