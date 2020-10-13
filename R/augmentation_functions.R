### sample_tree is a new version of the augment_tree function
### It can be used to simulate the extinct part of a tree, but also can be
### used to simulate a full tree
### The input is the object model which contains 


## sampling missing part of observed tree
sample_tree <- function(diversification_model,  # test1Ok
                        phylo,
                        max_num_species=100000){  
  
  # preparation

  tree = phylo2emph(phylo)
  ct = tree$ct
  cbt = 0
  
  while(cbt < ct){
    
    brts = c(tree$extant$brts,tree$extinct$brts,tree$ct)
    next_bt = min(brts[brts>cbt])

    # Draw speciation 
    next_speciation_time = draw_speciation(cbt,
                    next_bt,
                    diversification_model,
                    tree=tree,
                    full_tree=FALSE)
    if(next_speciation_time<next_bt){
      # resolve allocation
      allocation = draw_allocation(next_speciation_time,
                                 ct,
                                 diversification_model,
                                 tree)
 
      # update tree
      tree$extinct = rbind(tree$extinct,allocation)
        
      if(nrow(tree$extinct) > max_num_species){
          stop("Current parameters leds to a large number of species")
      }
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  return(tree)
}



draw_speciation <- function(cbt,
                            next_bt,
                            diversification_model,
                            tree=NULL,
                            full_tree=FALSE){
  
  if(full_tree){
    nsr = sum_of_rates
  }else{
    if(is.null(tree)){
      stop("The extant tree is needed as input")
    }
    nsr = nh_rate
  }
  key = 0 
  while(key == 0 & cbt < next_bt){
    
    lambda_max = optim(cbt,
                       fn = nsr,
                       tree = tree,
                       diversification_model = diversification_model,
                       lower = cbt,
                       upper = next_bt,
                       method ="L-BFGS-B",
                       control=list(fnscale=-1))$value
    
    u1 = runif(1)
    if(lambda_max==0){
      cbt = Inf
    }else{
      cbt = cbt - log(x = u1)/lambda_max
    }
    
    if(cbt < next_bt){
      u2 = runif(1)
      
      pt = nsr(tm=cbt,
               tree=tree,
               diversification_model=diversification_model)/lambda_max
      
      if(u2<pt){
        key = 1
      }
    }
  }
  
  return(cbt)
  
}

draw_allocation <- function(speciation_time,ct,diversification_model,tree){
  #draw extinction time
  extinction_time = speciation_time + 
        truncdist::rtrunc(1,
                          "exp",
                          a = 0,
                          b = (ct-speciation_time),
                          rate=diversification_model$pars[1])
  
  #choose the parent & child species
  
  current_species <- get_current_species(tm = speciation_time,
                                         tree = tree)
  parent <- sample(current_species,size = 1,
                   prob=speciation_rate(tm = speciation_time,
                                        tree = tree,
                                        diversification_model = diversification_model,
                                        sum_rates = FALSE))
  
  child <- max(c(tree$extant$child,tree$extinct$child))+1
  
  allocation = data.frame(brts=speciation_time,
                                 parent=parent,
                                 child=child,
                                 t_ext=extinction_time)
  
  return(allocation)
  
}

