## sample_tree is a new version of the augment_tree function
## It can be used to simulate the extinct part of a tree, but also can be
## used to simulate a full tree
## The input is the object model which contains 

sample_tree <- function(diversification_model,phylo=NULL,max_num_species=100000){  
  
  ## preparation

  tree = phylo2emph(phylo)
  ct = tree$ct
  cbt = 0
  
  while(cbt < ct){
    
    tree_extant = get_extant(tm,tree)
    brts = c(tree_extant$brts,ct)
    next_bt = min(brts[brts>cbt])

    ### Draw speciation 
    next_speciation_time = draw_speciation(cbt,
                    next_bt,
                    diversification_model,
                    tree=list(extant=tree_extant,
                              extinct=data.frame(brts=NULL,
                                                 parent=NULL,
                                                 child=NULL,
                                                 t_ext=NULL),
                              ct = ct),
                    full_tree=FALSE)
    if(next_speciation_time<next_bt){
      ## resolve allocation
      allocation = draw_allocation(next_speciation_time,
                                 ct,
                                 diversification_model,
                                 list(extant=tree_extant,
                                      extinct=data.frame(brts=NULL,
                                                         parent=NULL,
                                                         child=NULL,
                                                         t_ext=NULL),
                                      ct = ct))
 
      ## update tree

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
    nsr = nh_speciation_rate
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
      
      pt = nsr(cbt,
               tree,
               diversification_model)/lambda_max
      
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


draw_event_time <- function(cbt,
                            next_bt,
                            diversification_model,
                            tree){
  

  nsr = sum_of_rates

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
      
      pt = nsr(cbt,
               tree,
               diversification_model)/lambda_max
      
      if(u2<pt){
        key = 1
      }
    }
  }
  
  return(cbt)
  
}
#####

augment_tree <- function(phylo,
                         model){  
  pars = model$pars
  model = model$model 
  tree = phylo2emph(phylo)
  cbt = 0
  b = tree$ct
  mu = extinction_rate(tm = NULL,#constant extinction rate case
                       pars = pars,
                       model = model,
                       sum_rate = FALSE,
                       tree = NULL) 
  num_missing_branches=0
  
  while(cbt < b){
    brts = c(tree$extant$brts,tree$extinct$brts)
    next_bt = min(c(b,brts[brts>cbt]))
    
    lambda_max = optim(cbt,
                       fn = nh_rate,
                       pars = pars,
                       tree = tree,
                       model = model,
                       lower = cbt,
                       upper = next_bt,
                       method = "L-BFGS-B",
                       control = list(fnscale=-1)
    )$value
    
    if(lambda_max>500){
      stop("Current parameters leds to a huge speciation rate")
    }
    
    u1 = runif(1)
    next_speciation_time = cbt - log(x = u1)/lambda_max
    
    if(next_speciation_time < next_bt){  ## 
      u2 = runif(1)
      pt = nh_rate(x = next_speciation_time, 
                   model = model, 
                   pars = pars, 
                   tree = tree)/lambda_max
      if(u2 < pt){
        
        #draw extinction time
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        
        #choose the parent & child species
        
        current_species <- get_current_species(tm = next_speciation_time,
                                               tree = tree)
        parent <- sample(current_species,size = 1,
                         prob=speciation_rate(tm = next_speciation_time,
                                              tree = tree,
                                              pars = pars,
                                              model = model,
                                              sum_lambda = FALSE))
        
        child <- max(c(tree$extant$child,tree$extinct$child))+1
        
        
        # add to new tree
        to_add = data.frame(brts=next_speciation_time,
                            parent=parent,
                            child=child,
                            t_ext=extinction_time)
        tree$extinct = rbind(tree$extinct,to_add)
        
      }
      num_missing_branches <- num_missing_branches + 1
      if(num_missing_branches > 50000){
        stop("Current parameters leds to a large number of species")
      }
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  
  ## traits 
  tree$traits = list(n = n_for_all_bt(tree),
                     pd=2)
  ##
  
  return(tree)
}


augment_tree <- function(phylo,
                         model){  
  pars = model$pars
  model = model$model 
  tree = phylo2emph(phylo)
  cbt = 0
  b = tree$ct
  mu = extinction_rate(tm = NULL,#constant extinction rate case
                       pars = pars,
                       model = model,
                       sum_rate = FALSE,
                       tree = NULL) 
  num_missing_branches=0
  
  while(cbt < b){
    brts = c(tree$extant$brts,tree$extinct$brts)
    next_bt = min(c(b,brts[brts>cbt]))
    
    lambda_max = optim(cbt,
                       fn = nh_rate,
                       pars = pars,
                       tree = tree,
                       model = model,
                       lower = cbt,
                       upper = next_bt,
                       method = "L-BFGS-B",
                       control = list(fnscale=-1)
                       )$value
    
    if(lambda_max>500){
      stop("Current parameters leds to a huge speciation rate")
    }

    u1 = runif(1)
    next_speciation_time = cbt - log(x = u1)/lambda_max
    
    if(next_speciation_time < next_bt){  ## 
      u2 = runif(1)
      pt = nh_rate(x = next_speciation_time, 
                   model = model, 
                   pars = pars, 
                   tree = tree)/lambda_max
      if(u2 < pt){
        
        #draw extinction time
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        
        #choose the parent & child species
        
        current_species <- get_current_species(tm = next_speciation_time,
                                               tree = tree)
        parent <- sample(current_species,size = 1,
                         prob=speciation_rate(tm = next_speciation_time,
                                              tree = tree,
                                              pars = pars,
                                              model = model,
                                              sum_lambda = FALSE))
        
        child <- max(c(tree$extant$child,tree$extinct$child))+1
      
        
        # add to new tree
        to_add = data.frame(brts=next_speciation_time,
                            parent=parent,
                            child=child,
                            t_ext=extinction_time)
        tree$extinct = rbind(tree$extinct,to_add)
        
      }
      num_missing_branches <- num_missing_branches + 1
      if(num_missing_branches > 50000){
        stop("Current parameters leds to a large number of species")
      }
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  
  ## traits 
  tree$traits = list(n = n_for_all_bt(tree),
                     pd=2)
  ##

  return(tree)
}