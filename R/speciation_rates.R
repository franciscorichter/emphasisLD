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
  val = speciation_rate(tm = tm,
                  tree=tree,diversification_model = 
                  diversification_model,sum_rates = 
                  T) + 
        extinction_rate(tm,
                   tree,
                   diversification_model,
                   TRUE)
 
  return(val)
}

nh_rate <- function(tm,diversification_model,tree){
  val = speciation_rate(tm,
                  tree,
                  diversification_model,
                  sum_rates= TRUE)*
    (1-exp( -(tree$ct-tm) * mu.constant(tm,tree,diversification_model$pars,single_species = T)))
  return(val)
}

##########################################################
##### model-specific rates #################
##########################################################

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

lambda.ldpd <- function(tm,
                        tree,
                        pars,
                        sum_rates=FALSE){
  tree_extant = get_extant(tm,tree)
  gpd = GPD(tree_extant,tm)
  N = nrow(gpd)
  rates = pars[2] + pars[3]*N+ pars[4]*colSums(gpd)/(N-1)
  species = names(rates)
  lambdas = pmax(0,rates)
  names(lambdas) = species
  if(sum_rates) lambdas = sum(lambdas)
  return(lambdas)
}


lambda.ldpd2 <- function(tm,
                         tree,
                         pars,
                         sum_rates=FALSE){
  tree_extant = get_extant(tm,tree)
  gpd = GPD(tree_extant,tm)
  N = nrow(gpd)
  rates = pars[2] + pars[3]*N + (pars[4]*(N-1))/colSums(gpd)
  species = names(rates)
  lambdas = pmax(0,rates)
  names(lambdas) = species
  if(sum_rates) lambdas = sum(lambdas)
  return(lambdas)
}
# extinction rates

mu.constant <- function(tm,
                        tree,
                        pars,
                        sum_rates=FALSE,
                        single_species=FALSE){
  mu = max(0,pars[1])
  if(!single_species){
    N = sapply(tm, n_from_time,tree=tree)
    if(sum_rates){
      mu = N*mu
    }else{
      mu = rep(mu,N)
    }
  }
  return(mu)
}

