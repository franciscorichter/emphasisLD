############################################################

loglik.tree <- function(tree, diversification_model){
  to = tree$to
  to = head(to,-1)
  to[to!=0] = 1
  spec_times = tree$brts[c(to,0)==1]
  
  brts = c(tree$extant$brts,tree$extinct$brts)
  brts = brts[order(brts)]
  rates = sapply(brts, rate_at_bt, tree=tree, diversification_model=diversification_model)
  
  inte = intensity(tree, diversification_model)
  loglik = sum(log(rates)) + sum(log(extinctions)) - sum(inte)
  return(loglik)
}

rate_at_bt <- function(tm,tree,diversification_model){
  # oonly for rpd1 for the moment
  if(tm%in%tree$extinct$t_ext){
    val = extinction_rate(tm,tree,diversification_model)[1]
  }else{
    val = speciation_rate(tm,tree,diversification_model)[1]
  }
  return(val)
}

intensity <- function(tree, pars, model){
  nh_rate <- function(x){
    speciation_rate(tm=x,
                    tree = tree,
                    pars = pars,
                    model = model,
                    soc=tree$n[1],
                    sum_lambda = TRUE)+
      extinction_rate(tm=x,
                      tree = tree,
                      pars = pars,
                      model = model,
                      soc=tree$n[1],
                      sum_rate = TRUE)
  }
  brts_i = tree$brts
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),
                            xa = brts_im1[i]+0.00000000001,
                            xb = brts_i[i])
  }
  return(inte)
}
