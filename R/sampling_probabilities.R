# augmentaion (sampling) probability 

sampling_prob <- function(tree,pars,model){
  
  inte = intensity.numerical(tree = tree,pars = pars,model=model)

  if(nrow(tree$extinct)>0){
    top = head(tree$to,-1)
    N = tree$n
    soc = N[1]
    
    nb = N[missing_speciations]
    No = c(soc,soc+cumsum(top==2))[missing_speciations]
    Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
    
    brts_miss = tree$brts[missing_speciations]
    lambda_b = sapply(brts_miss,speciation_rate,tree = tree,pars = pars,model = model,soc=soc)
    text = tree$t_ext[missing_speciations]-brts_miss
    mu = max(0,pars[1])
    logg = -sum(inte)+sum(log(nb)+log(mu)+log(lambda_b)-mu*text-log(2*No+Ne))
  }else{
    logg = -sum(inte)
  }
  return(logg)
}


intensity.numerical <- function(tree, pars, model){
  brts = sort(c(tree$extant$brts,tree$extinct$brts,tree$extinct$t_ext,tree$ct))
  brts_i = brts[brts!=0]
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),
                            xa = brts_im1[i],
                            xb = brts_i[i],
                            model=model,
                            pars=pars,
                            tree=tree)
  }
  return(inte)
}


