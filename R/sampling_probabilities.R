# augmentaion (sampling) probability 

sampling_prob <- function(tree,pars,model){
  
  inte = intensity_missing(tree = tree,diversification_model)
  logg = -sum(inte)
  
  if(nrow(tree$extinct)>0){
    
    extant = tree$extant[tree$extant$brts!=0,]
    brts = c(extant$brts,tree$extinct$brts,tree$extinct$t_ext)
    to = c(rep(2,nrow(extant)),rep(1,nrow(tree$extinct)),rep(0,nrow(tree$extinct)))
    missing_speciations = (tree$to == 1)

    No = c(2,2+cumsum(to==2))[missing_speciations]
    Ne = c(0,cumsum(to==1)-cumsum(to==0))[missing_speciations]
    nb - No+Ne
    
    brts_miss = tree$brts[missing_speciations]
    sum_lambda_m = sapply(brts_miss,
                      speciation_rate,
                      tree = tree,
                      diversification_model=diversification_model,
                      sum_rates = TRUE)

    text = tree$extinct$t_ext-tree$extinct$brts
    mu = sapply(brts_miss, speciation_rate,
                tree = tree,
                diversification_model=diversification_model,
                sum_rates =FALSE)
    logg = logg+sum(log(mu)+log(sum_lambda_m)-mu*text-log(2*No+Ne))
  }

  return(logg)
}


intensity_missing <- function(tree, diversification_model){
  pars = diversification_model$pars
  model = diversification_model$model
  brts = sort(c(tree$extant$brts,tree$extinct$brts,tree$extinct$t_ext,tree$ct))
  brts_i = brts[brts!=0]
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),
                            xa = brts_im1[i],
                            xb = brts_i[i],
                            diversification_model,
                            tree=tree)
  }
  return(inte)
}


