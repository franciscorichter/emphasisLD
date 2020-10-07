############################################################

loglik.tree <- function(tree, diversification_model){
  
  brts = c(tree$extant$brts[-(1:2)],tree$extinct$brts,tree$extinct$t_ext)
  
  rates = sapply(brts, rate_at_bt, tree=tree, diversification_model=diversification_model)
  
  inte = intensity(tree, diversification_model)
  
  loglik = sum(log(rates)) + sum(log(extinctions)) - sum(inte)
  return(loglik)
}

rate_at_bt <- function(tm,tree,diversification_model){
  # only for rpd1 for the moment
  if(tm%in%tree$extinct$t_ext){
    val = extinction_rate(tm,tree,diversification_model)[1]
  }else{
    val = speciation_rate(tm,tree,diversification_model)[1]
  }
  return(val)
}

intensity <- function(tree, diversification_model){
  brts_i = sort(c(tree$extant$brts[-(1:2)],tree$extinct$brts,tree$extinct$t_ext))
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = sum_of_rates,
                            xa = brts_im1[i]+0.00000000001,
                            xb = brts_i[i],
                            tree = tree,
                            diversification_model = diversification_model)
  }
  return(inte)
}
